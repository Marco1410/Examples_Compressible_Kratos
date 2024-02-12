import numpy as np

####################################################
# CLASSES FOR DEFINING PARAMETRIZATION STRATEGIES
####################################################

# NOTE: parametrizations are to be plugged-in to one of the geometries defined in geom.py. 
class NACAparam:
  """Computes the geometry of NACA airfoil given c, p, m, t, or a 4 digit ID.
    Code adapted from https://people.math.sc.edu/Burkardt/py_src/naca/naca.html"""

  def __init__(self, m=0.0,p=0.0,t=0.0, NACAid=None):
    '''
    The NACA can be defined either by the three parameters (m,p,t) or by the "NACAid"  
    (nomenclature). If the NACAid is given, it is used to generate the geometry, otherwise,
    if it is None, the m,p,t paramters are employed.
    m, p, t: NACA values, given before normalization (ex. t=12.0 stands for 0.12 thicnkess)
    NACAid(string): 4 digits code of the NACA reference
    '''
    if NACAid==None:
      self.m = m / 100.0
      self.p = p / 10.0
      self.t = t / 100.0
    else:   
      self.NACAid=int(NACAid)
      # Compute airfoil characteristics
      self._get_airfoil_params_from_NACAid()

  def get_coordinates (self, x, c):
    '''
    Gives the geometry of NACA cambered 4-digit airfoil.
    Parameters:
       Input, real C, the chord length.
       Input, real X(*), points along the unit chord [0.0 <= X(*) <= 1.0].
       Output, real XU(*), YU(*), XL(*), YL(*), for each value of XC, measured
       along the camber line, the corresponding values (XU,YU) on the upper
       airfoil surface and (XL,YL) on the lower airfoil surface.
    '''
    xc=x*c

    i = np.nonzero ( 0.0 <= xc / c ) and np.nonzero ( xc / c <= self.p   )
    j = np.nonzero ( xc / c <= 1.0 ) and np.nonzero ( self.p   <= xc / c )
    k = np.nonzero ( xc / c < 0.0 ) or np.nonzero ( 1.0 < xc / c )

    n = len ( xc )

    divisor = np.zeros ( n )
    divisor[i] =         self.p  ** 2
    divisor[j] = ( 1.0 - self.p ) ** 2
    divisor[k] = 1.0

    dycdx = 2.0 * self.m * ( self.p - xc / c ) / divisor

    theta = np.arctan ( dycdx )
     
    yt = 5.0 * self.t * c * ( \
       0.2969 * np.sqrt ( xc / c ) \
       + (((( \
         - 0.1015 ) * ( xc / c ) \
         + 0.2843 ) * ( xc / c ) \
         - 0.3516 ) * ( xc / c ) \
         - 0.1260 ) * ( xc / c ) )

    yc = np.zeros ( n )
    yc[i] = self.m * ( xc[i]     ) * ( 2.0 * self.p - xc[i] / c       ) /         self.p   ** 2
    yc[j] = self.m * ( xc[j] - c ) * ( 2.0 * self.p - xc[j] / c - 1.0 ) / ( 1.0 - self.p ) ** 2
    yc[k] = 0.0

    xu = xc - yt * np.sin ( theta )
    yu = yc + yt * np.cos ( theta )
    xl = xc + yt * np.sin ( theta )
    yl = yc - yt * np.cos ( theta )

    return xu,yu,xl,yl

  def _get_airfoil_params_from_NACAid(self):
    '''
    Computes the characteristic values of a NACA airfoil based on the given code ([0,9999]).
      - Real M, the maximum camber, as a percentage of the chord length [0 <= M <= 1.0]
      - Real P, the relative distance of the occurrence of the maximum 
        camber from the beginning of the chord [0 <= P <= 1.0]
      - Real T, the maximum thickness relative to the chord length [0 <= T <= 1.0].
    '''
    code=self.NACAid
    if ( code < 0 or 9999 < code):
      print ( '' )
      print ( 'NACA4_MPT - Fatal error!' )
      print ( '  CODE should be an integer between 0 and 9999.' )
      exit ()

    m = ( code // 1000 )
    code = code - m * 1000
    self.m = m / 100.0

    p = ( code // 100 )
    code = code - p * 100
    self.p = p / 10.0

    self.t = code / 100.0