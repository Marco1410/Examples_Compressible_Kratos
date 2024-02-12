import numpy as np
from scipy.interpolate import pchip_interpolate
import os,meshio,gmsh,pathlib

print("INFO: working with meshio version from: %s" %(meshio.__file__))

####################################################
# CLASS FOR MESHING SIMPLE EXT. AERODYNAMIC CASES ##
####################################################

class extaero:
  """Manages the generation of 2D/2D.5 mesh files for external aero. A boundary layer around the datafile
  shape can be created. This shape will be built by performing an independent splining of the upper cloud 
  of points (xu,yu) and the lower cloud of points (xl,yl). 
  NOTE: depending on the geometry that one wants to generate, and the specified resolution, 
  building independent splines could break the continuity of the geometry"""

  def __init__(self, heightU=None, heightD=None, widthL=None, widthR=None, coneangle=0.0, xl=None, yl=None,
    xu=None, yu=None,roundmyTE=False,rotategeom=[]):
    '''Class initialization.
    * heightU, heightD, widthL, widthR: size of the CFD domain
    * If conangledeg [deg] is not 0.0, sets an opening angle in the direction of the mean flow
    * xl, yl, xu, yu: shape geometry (upper and lower sides)
    * roundmyTE: if set to True, it will assume the airfoil TE is provided as "open", and will
                create a round fillet to close it. Note that, if this flag is not active, the
                TE point of the upper and lower sided should be repeated 
                (see "_spline_top_bottom_lid" routine for full instructions).
    * rotategeom: is a list with [xrot,yrot,angledeg] to apply a rotation to xu,lu and xl,yl. 
                  Leave it as an empty list if no rotation is to be performed.
    '''
    self._set_width(widthL,widthR)
    if(coneangle!=0.0):
      angrad=np.deg2rad(coneangle)
      deltah=np.tan(angrad)*(widthL+widthR)
    else:
      deltah=0.0
    self._set_height(heightU,heightD,deltah)
    if(roundmyTE==True):
      # Creating a fillet and storing its width
      inixmax=max(max(xl),max(xu))
      xl,yl,xu,yu=_create_TE_fillet(xl,yl,xu,yu)
      endxmin=min(min(xl),min(xu))
      endxmax=max(max(xl),max(xu))
      ch=endxmax-endxmin
      TEl=endxmax-inixmax
      self.TEwidth=TEl/ch
    else:
      self.TEwidth=0.0
    self._set_datainfo(xl,yl,xu,yu)
    self.rotategeom=rotategeom
    # Sor far no mesh is there
    self.meshpath=None

  def generate_geom_and_mesh(self,outpath,lcwall=None,
                          lcedgemin=None,lcedgemax=None,edgedir=0,nblendspl=2,
                          blthick=None,blfirst=None,blratio=1.0,nfanmax=0,
                          threshold={},radialsizedict={},
                          extsize=None,mode="occ",algorithm=5,
                          planeid="xy",extrudel=None,extruden=None,
                          inletoutlet=False,topology="tri",
                          quadalgorithm="default",
                          isoptimize=False, meshfactor=1.0,
                          verbose=True, expertmode=True, outgeo=None):
    '''
    Generates a mesh in outpath. Extension is explicitely given by the user.
      * lcwall: Sets mesh size at the wall
      * lcedgemin: Same concept as lcwall, but this value will be used for the wall region
                 located at the minimum edge of the geometry (i.e. begining , e.g. LE). A blending
                 between lwall-lcedgemin will be produced, based on the direction 'edgedir'. As
                 the sizes of the splines can only be imposed at the end points, the code
                 will split the upper and lower splines into 'nblendspl', and keep imposing
                 the corresponding nodal sizes on the way.
                 Let it be None if no near-edge refinement is intended here.
      * lcedgemax: Same as lcedgemin, but deals with the max. edge (i.e. end, e.g. TE)
      * blthick: Boundary layer thickness. If set to None, no BL is generated.
                If you need additional granularity in size control, use also radialsizedict.
      * blfirst: Height of first boundary layer cell. If set to None, assumed to correspond to lcwall.
      * blratio: If blthick!=None, sets the expansion ratio of the boundary layer cells. If kept to 1.0
                all the cells are assumed to have the same size.
      * nfanmax: When doing boundary layers of e.g. airfoils, the TE may suffer from poor chordwise
                 density of grid lines. If that happens, use this integer to generate a singular
                 point around the TE and populate it with nfanmax lines.
      * threshold: is a dictionary with keys [distmin,distmax,lcmin], expressed as wall distance, that defines a background 
                 field for mesh size control. Below distmin, the size lcmin will be kept. Between distmin
                 and distmax, the size will transition from lcmin to extsize. Beyond that point, the
                 mesh size will be set to threshold["lcmin"]. Set distmin to 0.0 to start the size transition right from the wall.
                 If you do not want to use this feature, set it to {}.
                 extsize -                     /------------------
                                              /
                                             /
                                            /
                 threshold["lcmin"]  -o----/
                          |                |    |
                        Surface         distmin distmax
      * Use radialsizedict to define regions of refinement radially with respect to the mean position of the geometry. This allows to have a handle 
          on mesh size in the inner domain, so that e.g. additional regions of refinement can be set, complementing the BL refinement.
          Exemple: radialsizedict={1.0:0.5,2.0:7.0} will reinfornce the mesh size to be 0.5 at r=1.0, and 7.0 at r=2.0.
      * extsize: Set mesh size at the outer boundaries
      * mode: kernel of gmsh to use (geo -quite old format- or occ)
      * Choose algorithm according to http://gmsh.info/doc/texinfo/gmsh.html#index-Mesh_002eAlgorithm
          First test showed good performance for 5 -respects cell size without introducing too much anisotropy-
          and 8 -tries to make cells that can be combined later on to quads-. The gmsh default (i.e. 6) also give good results for most of the configurations.
      * Use planeid in case another reference frame is used in your solver. "xy" (default), and "xz" are supported.
      * extrudel: if =!None extrudes the mesh in the normal direction for the given value.
      * extruden: number of layers for the extrusion (only used if extrudel!=None)
      * inletoutlet: set it to True to identify inlet/outlet in the farfield, so that they are not combined as a single BC  
      * topology: sets the type of elements of the mesh, can get "tri", "quad" or "hybrid" 
                 (quads in the boundary layer -if any-, and triangles elsewhere). When doing 3D, these shapes are just extruded (prisms)
      * quadalgorithm: for some cases, such as boundary layers with very large blfirst, the standard recombination
              "default" algorithm (Blossom) may fail. For such cases, set quadalgorithm to "alternative". It leads to lower
              quality grids but can bypass Blossom failing. 
      * Set isoptimize to True to perform additional mesh optimizations (may have no effect)
      * meshfactor allows to scale all the imposed sizes on the mesh, allowing to generate coarser/finer grids
      * Set verbose==True to get info while generating the mesh.
      * If expertmode=True, then bypasses the WARNINGS targeting newbies (such as the fact that the 
          generated mesh will be too large, >10e4 nodes)
      * If outgeo is a string with a filename of a path (instead of None), dumps the geometry file there. Should have the extension .brep or .geo_unrolled.
    NOTE: Relies on direct interactions with gmsh-API for meshing [https://gitlab.onelab.info/gmsh/gmsh],
          insipired in example from examples/api/naca_boundary_layer_2d.py, and geo ideas from benchmarks/2d/naca12_2d.geo
    '''
    # Safety checks and gmsh initilization
    if(lcwall==None or extsize==None):
      print("ERROR: please provide a valid value for lcwall and extsize")
      exit()
    if(topology!="tri" and topology!="quad" and topology!="hybrid"):
      print("ERROR: Please set a valid topology for the mesh")
      exit()
    if(mode=="geo"):
      gmshgeomodel=gmsh.model.geo
    elif(mode=="occ"):
      gmshgeomodel=gmsh.model.occ
    else:
      print("ERROR: wrong mode chosen for mesh generation (options are geo or occ)")
      exit()
    if(planeid!="xy" and planeid!="xz"):
      print("ERROR: select a supported planeid (xy or xz)")
    self.gmshgeomodel=gmshgeomodel
    gmsh.initialize()

    # We save in the class the input variables
    self.outpath=outpath ; self.lcwall=lcwall  
    self.lcedgemin=lcedgemin ; self.lcedgemax=lcedgemax ; self.edgedir=edgedir ; self.nblendspl=nblendspl
    self.blthick=blthick ; self.blfirst=blfirst ; self.blratio= blratio ; self.nfanmax= nfanmax
    self.threshold=threshold ; self.radialsizedict=radialsizedict
    self.extsize=extsize ; self.mode=mode ; self.algorithm= algorithm
    self.planeid=planeid ; self.extrudel=extrudel ; self.extruden= extruden
    self.inletoutlet=inletoutlet ; self.topology= topology 
    self.quadalgorithm=quadalgorithm
    self.isoptimize=isoptimize ; self.meshfactor= meshfactor
    self.verbose=verbose ; self.expertmode=expertmode ; self.outgeo=outgeo

    # We track the total number of entities we are introducing 
    self.totnodes=0
    self.totcurves=0

    # We create our geometry by splining the upper and the lower surface, and introduce it
    self._compute_wall_geometry()
    self._rotate_walls_by_rotategeom()
    self._rotate_walls_toplane()
    self._introduce_wall_geometry()
    self._create_wall_curveloop()
    # Create outer domain
    self._compute_outer_geometry()
    self._rotate_outer_toplane()
    self._introduce_outer_geometry()    
    self._create_outer_curveloop()
    # Create planar surface between walls and outer domain
    self._create_plane_surface()
    # Adding thresold refinement if requested
    self._add_threshold()
    # Define boundary layer
    self._define_bl()
    # Adding radial refinement if requested.
    self._compute_radial_refinement()
    self._rotate_radial_refinement_toplane()
    self._introduce_radial_refinement()
    # Sets the BCs and extrude the geometry if needed
    self._set_bcs_and_extrude()
    # Configure gmsh options
    self._config_gmsh()

    # Write the .geo file if requested
    self._write_geo_file()
    # Generate the mesh
    self._generate_mesh()
    # Optimize mesh
    self.optimize_mesh()
    # Save mesh
    self._save_mesh()
    # Clear gmsh
    gmsh.clear()
    gmsh.finalize()

  def print_limits_x_info(self):
    '''Prints out the info of the outer CFD domain. x direction'''    
    print("INFO: CFD domain limits in direction 1: [%s,%s]" %(-self.widthL,self.widthR))

  def print_limits_y_info(self):
    '''Prints out the info of the outer CFD domain. y direction'''    
    print("INFO: CFD domain limits in direction 2: [%s,%s]" %(-self.heightD,self.heightU))

  def print_data_info(self):
    '''Prints out the info of the considered wall shape'''    
    print("INFO: wall limits:")
    print("\tIn direction 1: [%s,%s]" %(self.wxmin,self.wxmax))
    print("\tIn direction 2: [%s,%s]" %(self.wymin,self.wymax))     

  def print_mesh_info(self):
    '''Prints out the info of the last generated mesh'''
    if(self.meshpath==None):
      print("ERROR: generate a mesh before asking for its info")
      exit()
    meshr = meshio.read(self.meshpath)
    print("......... MESH REPORT BEGINS .........")
    print("Number of nodes: %s" %(meshr.points.shape[0])) 
    celldict=meshr.cells_dict
    for keyn in list(celldict.keys()):
      print("Number of %s cells: %s" %(keyn,celldict[keyn].shape[0]))
    print("......... MESH REPORT ENDS .........")
    return meshr.points.shape[0],celldict

  def _set_datainfo(self,xl,yl,xu,yu):
    '''Sets the cloud of points that will be used as walls.'''
    if((xl is None) or (yl is None)):
      return
    if((xu is None) or (yu is None)):
      return
    # Passing to 3D by assigning zeroes to the lower (datal)
    # and upper (datau) parts
    self.datal=np.zeros((len(xl),3))
    self.datal[:,0]=xl
    self.datal[:,1]=yl
    #
    self.datau=np.zeros((len(xu),3))
    self.datau[:,0]=xu
    self.datau[:,1]=yu
    #
    self._compute_data_limits()
    # Print the info of the geometry of the problem
    self.print_data_info()

  def _set_height(self,heightU,heightD,deltah):
    '''Sets height of the CFD domain. heightU for the top boundary,
    and heightD for the lower boundary. deltah will be used to stablish an opening angle'''
    self.heightUL=heightU-deltah*0.5
    self.heightUR=heightU+deltah*0.5
    self.heightDL=heightD-deltah*0.5
    self.heightDR=heightD+deltah*0.5
    self.heightU=max(self.heightUL,self.heightUR)
    self.heightD=max(self.heightDL,self.heightDR)
    self.print_limits_y_info()

  def _set_width(self,widthL,widthR):
    '''Sets width of the CFD domain. widthL for the left boundary, and
    widthR for the right boundary.'''
    self.widthL=widthL
    self.widthR=widthR
    self.print_limits_x_info()

  def _compute_data_limits(self):
    '''Computes the limits of the wall'''
    # Find limits of the wall
    self.wxmin = min(self.datal[:, 0].min(),self.datau[:, 0].min())
    self.wxmax = max(self.datal[:, 0].max(),self.datau[:, 0].max())
    self.wymin = min(self.datal[:, 1].min(),self.datau[:, 1].min())
    self.wymax = max(self.datal[:, 1].max(),self.datau[:, 1].max())

  def _compute_wall_geometry(self):
    '''
    Computes several splines (upper and lower), characterizing points and connectivity.
        self.datal: will be used as coordinates of the bootom lid
        self.datau: will be used as coordinates of the upper lid
    IMPORTANT rules
        Both curves should have a x2 common points (at the beggining and at the end)
            If a non-sharp TE is to be modeled, either include several artificial points around
            the TE or leave it open, and set "roundmyTE" to "True" when introducing the geometry
        Both curves should be defined in anti-clockwise direction
    A characteristic length to assign to the points (through self.lcwall, or a blending
        between self.lcdegemin/max and self.lcwall if requested).
    '''
    tol=1.0e-5
    # We do safety checks on the input
    for icoord in [0,1,2]:
      if(abs(self.datal[0,icoord]-self.datau[-1,icoord])>tol):
        print("ERROR:The initial vertex of the lower surface does not match")
        exit()
      if(abs(self.datau[0,icoord]-self.datal[-1,icoord])>tol):
        print("ERROR:The initial vertex of the upper surface does not match")
        exit()
    # Initialize what we will fill
    self.wallptsdict={} ; wallpointid=0
    self.wallsplsdict={} ; wallsplineid=0
    # Size of the initial data
    nl=len(self.datal[:,0])
    nu=len(self.datau[:,0])
    # We define the function that will help us to blend the sizes of the elements
    # so that additional refinement is found around the edges.
    # That will be a common blending function for upper/lower sides, using
    # as many control points as data points were provided
    nint=max(nl,nu)
    vmax=max(self.datal[:,self.edgedir].max(),self.datau[:,self.edgedir].max())
    vmin=min(self.datal[:,self.edgedir].min(),self.datau[:,self.edgedir].min())
    if(self.lcedgemin!=None or self.lcedgemax!=None):
      # Choosing the distribution we would like to have (moved to trapezoid to limit the refinement)
      # NOTE: Now not given as a variable, for simplicity. Also the trapezoidal margin is hardcoded
      distrib="trapez"
      #distrib="sine"
      if(distrib=="sine"):
        # distfact: distribution from 0 to 1 
        distfact=np.linspace (0.0, np.pi, nint)
        distfact=np.sin(distfact)
        # xc: corresponding x/c positions
        xc = np.linspace (vmin, vmax, nint)
      elif(distrib=="trapez"):
        # NOTE: Ensuring the transition happen in 5%, or forcing it to happen
        #       only within the fillet, if that was requested by the user.
        refsmooth=0.05
        deltatL=refsmooth
        print("INFO: Setting transition of lcedgemin/lcwall in %s percent chord" %(np.round(deltatL*100.0,2)))
        if(self.TEwidth==0.0):
          deltatR=refsmooth
          print("INFO: Setting transition of lcwall/lcdedgemax in %s percent chord" %(np.round(deltatR*100.0,2)))
        else:
          factsmooth=1.0
          deltatR=min(self.TEwidth*factsmooth,refsmooth)
          print("INFO: as TE was rounded, setting transition lcwall/lcdedgemax in %s percent chord" %(np.round(deltatR*100.0,2)))
        # Taking a given smoothclamp function 0-1 at both sides, and scale the abcissae based on input, from 0-1
        xcl,smoothL=_smooth_function(int(nint/3))
        xcl=xcl*deltatL
        xcr,smoothR=_smooth_function(int(nint/3),flipme=True)
        xcr=xcr*deltatR
        xcr=xcr+(1.0-deltatR)
        # Scaling with vmin/vmax
        xcr=xcr*(vmax-vmin)
        xcr=xcr+vmin
        # Fill center with remaining points (xcc does not dontain starting/ending point)
        Cnp=nint-len(smoothR)-len(smoothL)
        smoothC=np.array(Cnp*[1.0])
        distfact=np.concatenate((smoothL,smoothC,smoothR))
        xcc=np.linspace(xcl[-1],xcr[0],Cnp+1,endpoint=False)[1:] 
        xc=np.concatenate((xcl,xcc,xcr))
      # Scaling distfact with lcedgemin and lcedgemax
      # Split distfact in two halves
      midpoint = len(distfact) // 2
      bmin = distfact[:midpoint]
      bmax = distfact[midpoint:]
      # Build the final sizes array
      ycmin = self.lcedgemin-(self.lcedgemin-self.lcwall)*bmin
      ycmax = self.lcedgemax-(self.lcedgemax-self.lcwall)*bmax
      yc=np.concatenate((ycmin,ycmax))
    elif(self.lcedgemin==None and self.lcedgemax==None):
      # simply set lcedgemin and lcedgemax to lcwall
      xc = np.linspace (vmin, vmax, nint)
      yc = np.array(nint*[self.lcwall])
    # Deactivate blending if same value is found for lcedge and lcwall
    if(self.lcedgemin==self.lcwall and self.lcedgemax==self.lcwall):
      self.nblendspl=1
    # We define a structure with the geometry to walk (anti-clockwise)
    # IMPORTANT: As sizes along curves can only be properly imposed at the
    #            starting/ending points, we are forced to split the upper
    #            and lower surfaces into different splines (otherwise
    #            a double-sided smooting will never work)
    if((self.lcedgemin!=None or self.lcedgemax!=None) and self.nblendspl>1):
      spllst=[]
      for segment in _split_3d_data(self.datal,self.nblendspl):
        spllst.append({"datasp":segment})
      for segment in _split_3d_data(self.datau,self.nblendspl):
        spllst.append({"datasp":segment})
    else:
      spllst=[
      {"datasp" : self.datal},
      {"datasp" : self.datau},
      ]
    # We define a list with the points that we will introduce
    # IMPORTANT: We will NOT place two points at the same spot,
    #            so that we will start iterating from 1
    npspl=len(spllst)
    for curveid in range(0,npspl):
      # coordinates of this curve
      datasp=spllst[curveid]["datasp"]
      # Interpolate the sizes of the mesh
      resolnow=datasp[:,self.edgedir]
      sizes = pchip_interpolate(xc,yc,resolnow)
      # Save number of points in dictionary
      nlen=len(resolnow)
      spllst[curveid]["np"]=nlen
      for ipoint in range(1,nlen):
        pc=[datasp[ipoint,0],datasp[ipoint,1],datasp[ipoint,2]]
        lchar=sizes[ipoint]
        self.wallptsdict[wallpointid]={}
        self.wallptsdict[wallpointid]["coord"]=pc
        self.wallptsdict[wallpointid]["size"]=lchar
        wallpointid=wallpointid+1
    ptslst=list(range(wallpointid))
    # To introduce the splines, however, we will have to explicitely
    # refer to the common points. So that we rebuild the lists here
    npt=0
    for curveid in range(0,npspl):
      istart=npt
      nnow=spllst[curveid]["np"]
      iend=npt+nnow-1
      npt=iend
      if(curveid==0):
        # First curve takes initial point from last curve
        pts=[ptslst[-1]]+ptslst[istart:iend]
      else:
        pts=ptslst[istart-1:iend]
      spllst[curveid]["pts"]=pts
    # We are ready to loop and keep indtroducing the splines
    allsplines=[]
    for curveid in range(0,npspl):
      pts=spllst[curveid]["pts"]
      self.wallsplsdict[wallsplineid]=pts
      wallsplineid=wallsplineid+1

  def _rotate_walls_by_rotategeom(self):
    '''Rotates points and curves of the walls, using rotategeom (if requested)'''
    if(self.rotategeom==[]):
      return
    for pointid in list(self.wallptsdict.keys()):
      c=self.wallptsdict[pointid]["coord"]
      c=self._apply_rotategeom(c)
      self.wallptsdict[pointid]["coord"]=c

  def _rotate_walls_toplane(self):
    '''Rotates points and curves of walls to planeid (if requested)'''
    if(self.planeid=="xy"):
      return
    for pointid in list(self.wallptsdict.keys()):
      c=self.wallptsdict[pointid]["coord"]
      c=self._apply_rotateplaneid(c)
      self.wallptsdict[pointid]["coord"]=c

  def _introduce_wall_geometry(self):
    '''Introduces points of wall geometry, and creates splines based on connectivity'''
    # Initialize what we will fill
    self.wallpts=[] ; self.wallspls=[]
    iniptid=self.totnodes+1
    ptslst=list(self.wallptsdict.keys())
    for wallpointid in ptslst:
      pts=self.wallptsdict[wallpointid]
      ptid=wallpointid+iniptid
      size=pts["size"]
      c=pts["coord"]
      self.gmshgeomodel.addPoint(c[0],c[1],c[2],size,ptid)
      self.wallpts.append(ptid)
    #
    self.totnodes=self.totnodes+len(ptslst)
    self.gmshgeomodel.synchronize()
    #
    spllst=list(self.wallsplsdict.keys())
    inisplid=self.totcurves+1
    for wallsplineid in spllst:
      pts=self.wallsplsdict[wallsplineid]
      # shift by initial point id
      pts = [x + iniptid for x in pts]
      idspl=wallsplineid+inisplid
      # NOTE: addSpline can also work for some cases, but less general
      self.gmshgeomodel.addBSpline(pts,idspl)
      self.wallspls.append(idspl)
    self.totcurves=self.totcurves+len(spllst)

  def _create_wall_curveloop(self):
    '''Creates a curve loop from the wall curves'''
    cl=self.gmshgeomodel.addCurveLoop(self.wallspls)
    self.wallscloop=cl
    self.gmshgeomodel.synchronize()

  def _compute_outer_geometry(self):
    '''Compute points and lines of the outer domain'''
    # Initialize what we will fill
    self.outerptsdict={} ; self.outerlsdict={} 
    # Surface for numerical domain with a hole based on input geometry
    plst=[
      [-self.widthL, -self.heightDL, 0.0],
      [ self.widthR, -self.heightDR, 0.0],
      [ self.widthR,  self.heightUR, 0.0],
      [-self.widthL,  self.heightUL, 0.0]
    ]
    # introduce points
    for i in range(0,4):
      self.outerptsdict[i]={}
      self.outerptsdict[i]["coord"]=plst[i]
      self.outerptsdict[i]["size"]=self.extsize
    # introduce lines
    self.outerlsdict[0]=[0,1]
    self.outerlsdict[1]=[1,2]
    self.outerlsdict[2]=[2,3]
    self.outerlsdict[3]=[3,0]

  def _rotate_outer_toplane(self):
    '''Rotates points and curves of outer geometry to planeid (if requested)'''
    if(self.planeid=="xy"):
      return
    for pointid in list(self.outerptsdict.keys()):
      c=self.outerptsdict[pointid]["coord"]
      c=self._apply_rotateplaneid(c)
      self.outerptsdict[pointid]["coord"]=c

  def _introduce_outer_geometry(self):
    '''Introduces points of outer geometry, and creates lines based on connectivity'''
    # Initialize what we will fill
    self.outerpts=[] ; self.outersls=[] 
    iniptid=self.totnodes+1
    pointlst= list(self.outerptsdict.keys())
    for pointid in pointlst:
      pts=self.outerptsdict[pointid]
      ptid=pointid+iniptid
      size=pts["size"]
      c=pts["coord"]
      self.gmshgeomodel.addPoint(c[0],c[1],c[2],size,ptid)
      self.outerpts.append(ptid)
    self.totnodes=self.totnodes+len(pointlst)
    #
    self.gmshgeomodel.synchronize()
    #
    inicid=self.totcurves+1
    curveslt= list(self.outerlsdict.keys())
    for cid in curveslt:
      pts=self.outerlsdict[cid]
      # shift by initial point id
      pts = [x + iniptid for x in pts]
      idl=cid+inicid
      self.gmshgeomodel.addLine(pts[0],pts[1],idl)
      self.outersls.append(idl)
    self.totcurves=self.totcurves+len(curveslt)

  def _create_outer_curveloop(self):
    '''Creates the curve loop of the external boundary'''
    cl = self.gmshgeomodel.addCurveLoop(self.outersls)
    self.outercloop=cl
    self.gmshgeomodel.synchronize()

  def _create_plane_surface(self):
    '''Generates the planar surface between walls and the outer domain'''
    s = self.gmshgeomodel.addPlaneSurface([self.outercloop,self.wallscloop])
    self.wall2outersurf=s
    self.gmshgeomodel.synchronize()

  def _define_bl(self):        
    '''Defines a boundary layer is requested by the user'''
    self.gmshgeomodel.synchronize()
    if(self.blthick==None):
      return
    if(self.blfirst==None):
      print("INFO: size of first BL cell not speficied. Forced to be lcwall")
      self.blfirst=self.lcwall
    fieldn=len(gmsh.model.mesh.field.list())+1
    gmsh.model.mesh.field.add('BoundaryLayer',fieldn)
    gmsh.model.mesh.field.setNumbers(fieldn,'CurvesList',self.wallspls)
    gmsh.model.mesh.field.setNumber(fieldn,'Size',self.blfirst)
    gmsh.model.mesh.field.setNumber(fieldn,'Ratio',self.blratio)
    # 
    if(self.topology=="hybrid" or self.topology=="quad"):
      isquad=1
    else:
      isquad=0
    gmsh.model.mesh.field.setNumber(fieldn,'Quads',isquad)
    gmsh.model.mesh.field.setNumber(fieldn,'Thickness',self.blthick)
    # create a fan at the point furthest from origin (e.g. trailing edge), if requested
    if(self.nfanmax>0):
      fantptid=self._get_furthest_pt()
      gmsh.model.mesh.field.setNumbers(fieldn,'FanPointsList',[fantptid])
      gmsh.model.mesh.field.setNumbers(fieldn,'FanPointsSizesList',[self.nfanmax])
    gmsh.model.mesh.field.setAsBoundaryLayer(fieldn)

  def _get_furthest_pt(self):
    '''Returns the point identity of the point furthest from origin'''
    dlist=[]
    for ptid in list(self.wallptsdict.keys()):
      pcoord=self.wallptsdict[ptid]["coord"]
      # Distance to 0
      d=np.linalg.norm(pcoord)
      dlist.append(d)
    dlist=np.array(dlist)
    maxind = np.argmax(dlist)
    ptmaxid=self.wallpts[maxind]
    ptmaxcoord=self.wallptsdict[maxind]["coord"]
    print("INFO: Installing fan at point %s, with coordinates %s" %(ptmaxid,ptmaxcoord))
    return ptmaxid

  def _compute_radial_refinement(self):
    ''' Computes radial refinement if requested. Done by imposing the size in the points of discretized circles'''
    # Initialize radial points
    self.radialptsdict={} 
    irad=0
    for radnow in list(self.radialsizedict.keys()):
      sizenow=self.radialsizedict[radnow]
      # estimate the number of points to include in the circle, based on targeted size
      npts=int((2.0*np.pi*radnow)/sizenow)
      # give a bit of slack, as otherwise we will be over constraining the mesh (the points will be meshed)
      downsampl=5
      npts=int(float(npts)/float(downsampl))
      theta = np.linspace(0, 2.0*np.pi, npts, endpoint=False)
      # Calculate the x and y coordinates of a series of points where size will be reinforced
      # Centered in the mean position of the walls
      xcent=0.5*(np.mean(self.datau[:,0])+np.mean(self.datal[:,0]))
      ycent=0.5*(np.mean(self.datau[:,1])+np.mean(self.datal[:,1]))
      zcent=0.5*(np.mean(self.datau[:,2])+np.mean(self.datal[:,2]))
      xcyl=radnow*np.cos(theta)
      ycyl=radnow*np.sin(theta)
      zcyl=np.zeros(npts)
      cyl =np.array([xcyl, ycyl, zcyl]).T
      for i in range(0,npts):
        pcylx=cyl[i,0]+xcent
        pcyly=cyl[i,1]+ycent
        pcylz=cyl[i,2]+zcent
        self.radialptsdict[irad]={}
        self.radialptsdict[irad]["coord"]=[pcylx,pcyly,pcylz]
        self.radialptsdict[irad]["size"]=sizenow
        irad=irad+1

  def _rotate_radial_refinement_toplane(self):
    '''Rotates points and curves of radial refinement points to planeid (if requested)'''
    if(self.planeid=="xy"):
      return
    for pointid in list(self.radialptsdict.keys()):
      c=self.radialptsdict[pointid]["coord"]
      c=self._apply_rotateplaneid(c)
      self.radialptsdict[pointid]["coord"]=c
      
  def _introduce_radial_refinement(self):
    ''' Once the radial refinement is computed, embeds it in the surface'''
    self.radialpts=[]
    for pointid in list(self.radialptsdict.keys()):
      sizenow=self.radialptsdict[pointid]["size"]
      [pcylx,pcyly,pcylz]=self.radialptsdict[pointid]["coord"]
      self.gmshgeomodel.synchronize()
      ps1=self.gmshgeomodel.addPoint(pcylx,pcyly,pcylz,sizenow)
      self.radialpts.append(ps1)
      self.gmshgeomodel.synchronize()
      itemdim=0 ; hostdim=2 
      gmsh.model.mesh.embed(itemdim,[ps1],hostdim,self.wall2outersurf)

  def _add_threshold(self):
    ''' Adds thresold if requested. Wall distance -based refinement based 
    on t10.py example of gmsh'''
    if(self.threshold=={}):
      return
    self.gmshgeomodel.synchronize()
    fieldn=len(gmsh.model.mesh.field.list())+1
    gmsh.model.mesh.field.add("Distance", fieldn)
    gmsh.model.mesh.field.setNumbers(fieldn, "PointsList", self.wallpts)
    gmsh.model.mesh.field.setNumbers(fieldn, "CurvesList", self.wallspls)
    # nsamp: number of points that will be used to check the distance field
    nsamp=500
    gmsh.model.mesh.field.setNumber(fieldn, "Sampling", nsamp)
    gmsh.model.mesh.field.add("Threshold", fieldn+1)
    gmsh.model.mesh.field.setNumber(fieldn+1, "InField", fieldn)
    # NOTE: hardcoded min/max values, as at the end of the day we will use the "min" operator
    gmsh.model.mesh.field.setNumber(fieldn+1, "SizeMin", self.threshold["lcmin"])
    gmsh.model.mesh.field.setNumber(fieldn+1, "SizeMax", self.extsize)
    gmsh.model.mesh.field.setNumber(fieldn+1, "DistMin", self.threshold["distmin"])
    gmsh.model.mesh.field.setNumber(fieldn+1, "DistMax", self.threshold["distmax"])
    gmsh.model.mesh.field.setAsBackgroundMesh(fieldn+1)

  def _set_bcs_and_extrude(self):
    '''Sets the BCs and extrude the geometry if needed'''
    self.gmshgeomodel.synchronize()
    if(self.extrudel!=None):
      # NOTE: we decompose the outcome of the extrusion based on the knowledge of the command.
      #       it is given by tuples, so that we access the second position.
      # NOTE: recombine is used, to create elements with edges parallels to extrusion,
      #       but this could be eventually be passed as a flag
      self.bottomsurf=self.wall2outersurf
      extruded=self._extrude_entity(2,[self.wall2outersurf])
      self.topsurf=extruded[0][1]
      self.volume=extruded[1][1]
      # NOTE: this decomposition assumes that we have a x4 sides domain.
      #       Additional coding is needed if more outer shapes are to be implemented
      self.p0surf=extruded[2][1]
      self.p1surf=extruded[3][1]
      self.p2surf=extruded[4][1]
      self.p3surf=extruded[5][1]
      self.wallsurf=[]
      for i in range(6,len(extruded)):
        self.wallsurf.append(extruded[i][1])
      #
      self._add_gmsh_bcs([self.bottomsurf],dim=2,name="bottomsurf")
      self._add_gmsh_bcs([self.topsurf],dim=2,name="topsurf")
      #
      self._add_gmsh_bcs([self.volume],dim=3,name="fluid")
      # Apply inlet/outlet if needed, and specify the rest of the farfield
      if(self.inletoutlet==True):
        self._add_gmsh_bcs([self.p0surf],dim=2,name="farfielddown")
        self._add_gmsh_bcs([self.p1surf],dim=2,name="outlet")
        self._add_gmsh_bcs([self.p2surf],dim=2,name="farfieldup")
        self._add_gmsh_bcs([self.p3surf],dim=2,name="inlet")
      else:
        polysurf=[self.p0surf]+[self.p1surf]+[self.p2surf]+[self.p3surf]
        self._add_gmsh_bcs(polysurf,dim=2,name="farfield")
      #
      self._add_gmsh_bcs(self.wallsurf,dim=2,name="walls")
      #
    else:
      # Set BCs
      self._add_gmsh_bcs(self.wallspls,dim=1,name="walls")
      if(self.inletoutlet==True):
        # NOTE: Will only work if the outer domain has 4 sides (hardcoded indices)
        self._add_gmsh_bcs([self.outersls[1]] ,dim=1,name="outlet")
        self._add_gmsh_bcs([self.outersls[3]] ,dim=1,name="inlet")
        self._add_gmsh_bcs([self.outersls[0]]+[self.outersls[2]],dim=1,name="farfield")
      else:
        self._add_gmsh_bcs(self.outersls,dim=1,name="farfield")
      self._add_gmsh_bcs([self.wall2outersurf],dim=2,name="fluid")

  def _extrude_entity(self,ndim,entitylst,recombine=True):
    '''Extrudes the input entity, with dimension ndim, and returns the outcome.
       By default, recombines the elements in the direction of the extrusion'''
    if(self.planeid=="xy"):
      nvec=[0.0,0.0,self.extrudel]
    elif(self.planeid=="xz"):
      nvec=[0.0,self.extrudel,0.0]
    tuplelst=_get_tupledim(ndim,entitylst)
    # NOTE: one could also use "heights" to control the mesh size
    #       of the elements while extruding.
    outent=self.gmshgeomodel.extrude(tuplelst,nvec[0],nvec[1],nvec[2],
            numElements=[self.extruden],recombine=recombine)
    return outent

  def _add_gmsh_bcs(self,bcs,dim,name):
    '''Add physical group to gmsh, after synchonizing'''
    self.gmshgeomodel.synchronize()
    gmsh.model.addPhysicalGroup(dim=dim,tags=bcs,tag=-1,name=name)

  def _config_gmsh(self):
    '''Configures gmsh for our problem'''
    # IMPORTANT: avoid gmsh collapsing the quads into triangles when extruding
    gmsh.option.setNumber("Geometry.Tolerance",1e-15)
    # Setting all quads if requested
    if(self.topology=="quad"):
      gmsh.option.setNumber("Mesh.RecombineAll", 1)
      # Changing combination algorithm if requested
      if(self.quadalgorithm=="alternative"):
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm",0)
    # Overall control of mesh size
    minsizelst=[self.blfirst,self.lcwall,self.lcedgemin,self.lcedgemax] 
    minsize=min(x for x in minsizelst if x is not None)
    print("INFO: Allowing the mesh size to vary between %s and %s" %(minsize,self.extsize))
    gmsh.option.setNumber("Mesh.MeshSizeMin", minsize)
    gmsh.option.setNumber("Mesh.MeshSizeMax", self.extsize)

    # Including by-default smoothing of the mesh
    gmsh.option.setNumber("Mesh.Smoothing", 3)

    # Trying to ensure the repeatibility of the mesh
    # IMPORTANT: It did not work with current gmsh versions at least, as
    #            there seems to be some time-stamp dependent random generatorn
    #            that cannot be forced from current API.
    gmsh.option.setNumber("Mesh.RandomFactor", 1e-9)
    gmsh.option.setNumber("Mesh.RandomFactor3D", 1e-12)
    gmsh.option.setNumber("Mesh.RandomSeed", 1.0)

    # ...... Fine Tuning gmsh execution .....
    #   NOTE: we enable multithread for mesh creation/ for boolean operation
    #   IMPORTANT: Does not seem to work in all setups [depends on gmsh installation apparently]
    gmsh.option.setNumber("General.NumThreads", 8)
    gmsh.option.setNumber("Mesh.MaxNumThreads1D", 8)
    gmsh.option.setNumber("Mesh.MaxNumThreads2D", 8)
    gmsh.option.setNumber("Mesh.MaxNumThreads3D", 8)
    gmsh.option.setNumber("Geometry.OCCParallel", 1)
    
    # Forcing to check points size (as used in radialsizedict and while building surface splines)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 1)
    # Activate expert mode to bypass WARNINGS, if requested
    if(self.expertmode==True):
      gmsh.option.setNumber("General.ExpertMode",1)
       
    # Direct handles on gmsh setup (use for debug)
    # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # gmsh.option.setNumber("Mesh.MeshSizeFromCurvatureIsotropic", 0)
    # gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    # gmsh.option.setNumber("Mesh.MeshSizeFromParametricPoints", 1)
    # gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    if(self.meshfactor!=1.0):
      gmsh.option.setNumber("Mesh.MeshSizeFactor", self.meshfactor)

  def _write_geo_file(self):
    '''Writes a raw version of the geo file, for debugging. It could be eventually read by gmsh'''
    if(self.outgeo==None):
      return
    pathlst=pathlib.Path(self.outgeo).suffixes[-1]
    if((".brep" not in pathlst) and (".geo_unrolled" not in pathlst)):
      print("ERROR: only .brep and .geo_unrolled extensions allowed when dumping mesh geometry")
      exit()
    self.gmshgeomodel.synchronize()
    gmsh.write(self.outgeo)

  def _generate_mesh(self):
    '''Generates the mesh, once the problem is properly configured'''
    self.gmshgeomodel.synchronize()
    print("INFO: Generating mesh (non-deterministic process)")
    # Setting verbosity
    # 0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5: +status, 99: +debug
    if(self.verbose==False):
      verbosity=0
    elif(self.verbose==True):
      verbosity=5
    # Set dimensions of the problem
    if(self.extrudel!=None):
      ndim=3
    else:
      ndim=2
    for algotab in ["Mesh.Algorithm","Mesh.Algorithm3D"]:
      gmsh.option.setNumber(algotab, self.algorithm)
    gmsh.option.setNumber("General.Verbosity", verbosity)
    gmsh.model.mesh.generate(ndim)

  def optimize_mesh(self):
    '''Optimizes the mesh if requested. May have no effect'''
    if(self.isoptimize!=True):
      return
    # NOTE: this is the default option, that aims for tetrahedrals 
    algo=''
    if(self.extrudel==None):
      algo="Relocate2D"
    else:
      algo="Relocate3D"
    if(self.topology=="quad"):
      algo="QuadQuasiStructured"
    gmsh.model.mesh.optimize(algo,force=True)
  
  def _save_mesh(self):
    '''Saves already generated mesh'''
    print("INFO: Saving mesh") 
    # Save the mesh with the requested extension and record the path
    gmsh.write(self.outpath)
    self.meshpath=self.outpath

  def _apply_rotategeom(self,c):
    '''Apply requested geometrical rotation to wall-related points'''
    # NOTE: This rotation is done manually since gmsh API was found to be shady (does
    #        not keep entities' integrity)
    if(self.rotategeom==[]):
      return
    xrot=self.rotategeom[0]
    yrot=self.rotategeom[1]
    rotatedeg=self.rotategeom[2]
    # NOTE: negative so a that a positive rotategeom leads to pitch up
    angrad=-np.deg2rad(rotatedeg)
    Mz=np.matrix([[ np.cos(angrad), -np.sin(angrad), 0.0],
                  [ np.sin(angrad), np.cos(angrad), 0.0],
                  [ 0.0, 0.0, 1.0]])
    cout=_rotate_l(c,Mz,pref=[xrot,yrot,0.0])
    return cout

  def _apply_rotateplaneid(self,c):
    '''Applies a rotation to set the point to a non-standard plane'''
    # NOTE: This rotation is done manually since gmsh API was found to be shady (does
    #        not keep entities' integrity)
    if(self.planeid=="xy"):
      cout=c
    elif(self.planeid=="xz"):
      angrad=np.deg2rad(90.0)
      Mx=np.matrix([[ 1.0, 0.0          , 0.0         ],
                   [ 0.0, np.cos(angrad),-np.sin(angrad)],
                   [ 0.0, np.sin(angrad), np.cos(angrad)]])
      cout=_rotate_l(c,Mx)
    return cout

##############################################
##### METHODS NOT IN CLASS   #################
##############################################

def invert_mdpa_lines(infile,outfile,intag,outtag):
  '''
  NOTE: This method allows to remove the Kratos WARNING about the orientation
  of the conditions. For the particular airfoil examples included in this repo,
  that corresponds to a call like:
      cfdd.mesh.invert_mdpa_lines(kratospathinverted,kratospathfixed,
          "GUI group identifier: walls","End Conditions")
  However, this correction does not seem to have any effect on the results, and even
  some GiD-prepared models do have the same issue. So it is better to spare it.
  '''
  fin=open(infile,'r')
  fout=open(outfile,'w')
  inblock=False
  line = fin.readline()
  fout.write(line)
  while line:
    line = fin.readline()
    if((inblock==True) and (outtag in line)):
      inblock=False
    if(inblock==True):
      linespl=line.split()
      fout.write("%s %s %s %s\n" %(linespl[0],linespl[1],
                                  linespl[3],linespl[2]))
    else:
      fout.write(line)
    if((inblock==False) and (intag in line)):
      inblock=True
  fin.close()
  fout.close()

def from_gmsh41_to_gmsh22(outpath,binary=False):
  '''Passess from gmsh 4.1 to 2.1 format (ASCII)'''
  print("INFO: Pass from Gmsh4.1 to Gmsh2.2 for: %s" %(outpath))
  mesh = meshio.read(outpath)
  meshio.write(outpath,mesh,file_format="gmsh22",binary=binary)

def from_gmsh41_to_Kratosmdpa(inpath,outpath,removesource=True,
  dictbcs={"walls":"LineCondition2D2N","farfield":"LineCondition2D2N","fluid":"Element2D3N"}):
  '''Passess from gmsh 4.1 (inpath) to Kratos mdpa format (outpath). The dictionary
  dictbcs is used to re-name the conditions and elements in the mdpa, from geometry-based to
  formulation-based'''
  print("INFO: Pass from Gmsh4.1 to Kratos/mdpa for: %s" %(inpath))
  mesh = meshio.read(inpath)
  meshio.write(outpath,mesh,file_format="mdpa")
  # Remove original file if requested by the user
  if(removesource==True):
    os.remove(inpath)
  # Adapt conditions to what Kratos expects
  _adaptnamingkratosmdpa(outpath,dictbcs)

def _adaptnamingkratosmdpa(filepath,dictbcs):
  '''Renames the mdpa sections generated by meshio. This process is quite manual,
  and only takes into account the values passed in the dict. Should be exchanged by something else
  as the meshio mdpa interface (master) becomes updated'''
  fin=open(filepath,"r")
  data=fin.readlines()
  fin.close()
  fout=open(filepath,"w")
  dictkeys=list(dictbcs.keys())
  for line in data:
    linelst=line.rstrip('\n').split()
    linelst=[s.strip() for s in linelst]
    # We just check there is a common element, and remove the
    # 3rd element accoding to the input dict.
    # We also skip the last definitions, where the lines only have three letters
    # NOT BULLET PROOF!-> Just a temporary workaround
    if(len(linelst)>3  and _common_member(linelst,dictkeys)):
      lastword=linelst[-1]
      linelst[2]=dictbcs[lastword]
    linestr=' '.join(linelst)
    fout.write("%s\n" %(linestr))
  fout.close()

def _common_member(a, b):
  '''Checks if two lists have a common member'''
  if len(list(set(a) & set(b)))>0:
    return True
  return False

def _split_3d_data(data,ns):
  '''
  Splits data into ns segments, assuming a point is shared 
  by each of them
  '''
  out=[]
  ntot=len(data)
  nstep=int(ntot/ns)
  for i in range(0,ns):
    istart=i*nstep
    iend=min((i+1)*nstep+1,ntot+1)
    out.append(data[istart:iend,:])
  return out

def _create_TE_fillet(xl,yl,xu,yu):
  '''Assumes the input defines the lower and upper surfaces of a e.g. open TE airfoil,
     and creates a round fillet between them'''
  # find the center and radius
  p1=np.array([xu[0],yu[0]])
  p2=np.array([xl[-1],yl[-1]])
  c=0.5*(p1+p2)
  d=np.linalg.norm(p2-p1)
  r=0.5*d
  # find the last point
  u=(p1-p2)/d
  un=np.array([u[1],-u[0]])
  pf=c+un*r
  # Define our angular step [deg]
  dalpha=1.0
  # fill upper
  ang=np.arange(0.0, 90.0, dalpha)
  ang=np.deg2rad(ang)
  xpts=c[0]+r*np.cos(ang)
  ypts=c[1]+r*np.sin(ang)
  xu=np.concatenate((xpts,xu))
  yu=np.concatenate((ypts,yu))
  # fill lower
  ang=np.arange(0.0, 90.0, dalpha)
  ang=np.deg2rad(-ang)
  ang=np.flip(ang)
  xpts=c[0]+r*np.cos(ang)
  ypts=c[1]+r*np.sin(ang)
  xl=np.concatenate((xl,xpts))
  yl=np.concatenate((yl,ypts))
  #
  return xl,yl,xu,yu

def _smooth_function(npoints,flipme=False):
  '''Creates a smooth function from 0 to 1, both in x and y, using npoints. Flips the result is needed'''
  x = np.linspace(0, 1, npoints)
  # smooth version
  #y = -2*x**3 + 3*x**2
  # getting the start more flat
  y = -4*x**5 + 10*x**3
  if(flipme==True):
    y=np.flip(y)
  return x,y

def _get_tupledim(val,listin):
  '''Generate a list of tuples (val,item), where item is each of the items
  of the input list'''
  lst=[]
  for itemn in listin:
    lst.append((val,itemn))
  return lst 

def _rotate_l(cin,M,pref=[0.0,0.0,0.0]):
  '''Performs a rotation of a 3d point, given as a list. A ref. point for the rotation can be specified'''
  # shift by reference
  crot=[cin[0]-pref[0],cin[1]-pref[1],cin[2]-pref[2]]
  # perform rotation
  pta=np.dot(M,crot)
  pta=pta.tolist()[0]
  # undo shifting
  cout=[pta[0]+pref[0],pta[1]+pref[1],pta[2]+pref[2]]
  return cout


