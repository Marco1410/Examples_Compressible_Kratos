import numpy as np
from scipy.interpolate import interp1d,Akima1DInterpolator

####################################################
# CLASSES FOR GENERATING SIMPLE GEOMETRIES #########
####################################################

# NOTE: The types of geometry shown here require a parametrization as an input (to be included in param.py),
#       accounting for the method get_coordinates.
class airfoil:
    """Handles the generation and manipulation of airfoil geometries. Supports several parametrizations"""  
    def __init__(self, param=None, TEtype="base"):
      '''
      param: an airfoil parametrization class
      TEtype: type of LE to generate
            "base": generates PS and SS based on originally given data, without repeating any point
            "blunt": creates an additional point, repeated for SS and PS, at the mean position of
                     the corresponding original TE positions
            "sharp": gives a sharp TE, which position is estimated numerically.
      '''
      if param==None:
        print("ERROR: give a parametrization class to start an airfoil instance")
        exit()
      self.param=param
      self.TEtype=TEtype
    def generate_geometry(self, n, c=1.0, distri="sine", redist=True, rechord=False):
      ''' 
      Generates the PS and SS geometries of the requested airfoil.
        - n is an integer with the number of points that we want to use as an input
          for the parametrization of the surfaces (PS and SS). NOTE: is redist!=None,
          the output will not correspond to these values, even if TEtype=="base".
        - c is a float (chord) to scale the geometry
        - distri can take "equi" (equidistant in chord axis) or "sine" (will
          stack more point close to LE and TE)
      redist: this option allows to re-cluster the points based on requested distribution [distri], so that:
            False: the computed x coordinates by the parametrization will be returned, 
                   (adding a trailing edge point if needed), but no redistribution will be performed
            True: re-produces the original distribution [distri] on the new geometry.
      rechord: if set to True, this option rescales the airfoil so that the x coordinate goes
               from 0 to c. It will have an effect as long as a TEtype rather than "base" has been set,
               or if it is a non-symmetric airfoil.
      '''
      # Get requested distribution 
      x=get_axis_distribution(distri,n) 
      # Once we have the coordinate, we evaluate our airfoil
      xu,yu,xl,yl=self.param.get_coordinates(x,c)
      # We add a point for blunt or sharp TE if requested
      islong=False
      if(self.TEtype=="blunt"):
        islong=True
        yTE=0.5*(yu[-1]+yl[-1])
        # NOTE: One could also add a user-defined distance for XTE (instead of
        #       1e-4), but that falls back into "sharp"
        xTE=1.0001*c
      elif(self.TEtype=="sharp"):
        islong=True
        # Create splines of the surfaces, allowing for
        # extrapolation to be able to find their intersection
        # NOTE: we use half of the points to speed up the process
        inii=int(n/2)
        splu = interp1d(xu[inii:-1],yu[inii:-1],
          fill_value="extrapolate",kind="cubic")
        spll = interp1d(xl[inii:-1],yl[inii:-1],
          fill_value="extrapolate",kind="cubic")
        # Find intersection at the vicinity of LE
        xTE,yTE=find_intersection_splines(splu,spll,0.8*c,1.1*c)
        if(len(xTE)!=1):
          print("ERROR: unable to find sharp trailing edge location numerically")
          exit()
        xTE=xTE[0] ; yTE=yTE[0]
      elif(self.TEtype!="base"):
        print("ERROR: TEtype value not accepted when creating geometry")
        exit()
      # We add the TE point if requested
      if(islong):
        xl=np.append(xl,xTE) ; yl=np.append(yl,yTE)  
        xu=np.append(xu,xTE) ; yu=np.append(yu,yTE)
      # We re-check if a sharp TE airfoil was generated (as it could be that the "base" parametrization
      # already produced a sharp one)
      isharp=False
      if(xl[-1]==xu[-1] and yl[-1]==yu[-1]):
        isharp=True
      # We redistribute geometry if requested
      if(redist):
        # NOTE: for non-symmetric airfoils, the curvature of the LE may impose
        #       different trends on SS or PS (increasing-decreasing), that breaks
        #       the spline interpolation. We correct this by re-generating the curves
        if(np.min(xu)<0.0 or np.min(xl)<0.0):
          # Concatenate and avoid duplication of LE point
          xt=np.append(np.flip(xl,0),xu[1:])
          yt=np.append(np.flip(yl,0),yu[1:])
          ind=np.where(xt==xt.min())
          if(len(ind)!=1):
            print("ERROR: cannot re-build SS and PS in non-symmetric airfoil")
            exit()
          else:
            ind=int(ind[0])
          xl=xt[0:ind+1] ; yl=yt[0:ind+1]
          xu=xt[ind:]    ; yu=yt[ind:]
          xl=np.flip(xl) ; yl=np.flip(yl)
        # NOTE: Akima splines will be used to avoid near LE overshooting.
        #       These were not used while finding TE position as they do not allow
        #       for extrapolation.
        splu = Akima1DInterpolator(xu,yu)
        spll = Akima1DInterpolator(xl,yl)
        # NOTE: distributions will be different with open TE
        xu=xu[0]+x*(xu[-1]-xu[0]) 
        xl=xl[0]+x*(xl[-1]-xl[0]) 
        yu=splu(xu) ; yl=spll(xl)
        # Ensuring same mathematical point at TE when forcing closure
        # (precicion-based correction)
        if(isharp):
          yTEm=0.5*(yu[-1]+yl[-1])
          yu[-1]=yTEm ; yl[-1]=yTEm
        #
      if(rechord):
        # NOTE: just taking the max. when working with open TE. 
        if(isharp):
        	d=xu[-1]-xu[0]
        else:
        	d=max(xu[-1],xl[-1])-xu[0]
        fc=c/d
        xu=xu-xu[0] ; xl=xl-xl[0]
        xu=fc*xu ; yu=fc*yu
        xl=fc*xl ; yl=fc*yl
      # We already flip the upper surface to 
      # ensure the continuity
      self.xu=np.flip(xu,0) ; self.xl=xl 
      self.yu=np.flip(yu,0) ; self.yl=yl
      # Assemble them concatenated without LE duplication
      self.xa=np.append(self.xu,self.xl[1:])
      self.ya=np.append(self.yu,self.yl[1:])

    def write_airfoil(self,outfile,surftype="all"):
      '''Writes coordinates of the airfoil to a file.
         surftype: 
          all: both PS and SS surfaces
          PS: pressure side
          SS: suction side
       '''
      x,y=self.give_coordinates(surftype)
      np.savetxt(outfile, np.c_[x,y])   

    def give_coordinates(self,surftype="all"):
      '''Gets coordinates of the airfoil.
         surftype: 
          all: both PS and SS surfaces
          PS: pressure side
          SS: suction side
      '''
      if(surftype=="all"):
        x=self.xa ; y=self.ya
      elif(surftype=="SS"):
        x=self.xu ; y=self.yu
      elif(surftype=="PS"):
        x=self.xl ; y=self.yl
      else:
        print("ERROR: for airfoil coordinate requests, enter -all-, -PS- or -SS-")
        exit()
      return x,y

########################################
######## COMMON METHODS ################
########################################

def get_axis_distribution(distri,n):
  '''
  Distributes n points along an axis, from 0 to 1.0. returns result as array
  Supports equidistant and sine distributions (dirstri variable)
  '''
  if(distri!="equi" and distri!="sine"):
    print("ERROR: please enter a valid distri to generate airfoil geometry (i.e equi or sine)")
    exit()
  if(distri=="equi"):
    x=np.array(n*[1.0])
  elif(distri=="sine"):
    x=np.linspace (0.0, np.pi, n)
    x=np.sin(x)
    # Note: we extrapolate the very last point based on the last value
    #       otherwise we will place two points at the same spot
    x[-1]=x[-2]
  # We accumulate the segments to get the x coordinate and scle them to 1.0
  x=np.cumsum(x)
  x=x-x[0]
  x=x/x[-1]
  return x

def find_intersection_splines(spl1,spl2,xmin,xmax,nrec=1000):
    '''
    Search for the intersection of splines spl1 and spl2 in range [xmin,xmax]
    Uses nrec to create an array in the given limits
    Note: the spline should be defined in this range (i.e. extrapolation might be needed)
    '''
    #
    xeval=np.linspace(xmin,xmax, num=nrec)
    f=spl1(xeval)
    g=spl2(xeval)
    #
    idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
    return xeval[idx], f[idx]
