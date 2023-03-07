import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline as RBS
import matplotlib.pyplot as plt

class InertialLiftForceHelperVST(object):
    """This class helps with parsing/utilising computed inertial lift  
    force and migration data for a neutrally buoyant spherical  
    particle suspended in flow through curved ducts having a vertically
    symmetric trapezoidal cross-section and aspect ratio 2:1 or 4:1
    (measured as the width compared to the average/centre height).
    
    The class will interpolate the raw data for a desired bend radius of
    the duct. Interpolation of particle size is not currently supported.
    As such one must pick a particle size matching the provided data.
        
    Note: use of this helper class requires the files 
    'VST_data/strap2x1_lift_data.npz'
    'VST_data/strap4x1_lift_data.npz'
    to be located relative to the working directory.
    
    Please cite our research if you use this code/data.
    The model is developed in a JFM paper (doi.org/10.1017/jfm.2019.323).
    The dynamics of trapezoidal ducts is explored in (TODO: add details).
    
    This code is provided under an MIT license (see https://opensource.org/licenses/MIT).
    However, I ask that the associated research data not be directly distributed to
    third parties, rather I would appreciate if it was cloned directly from my respository.
    Please don't hesitate to contact me if you have any questions/queries.
    
    Brendan Harding, last updated 2023."""
    def __init__(self,a,R,ar=2,D=0):
        """Initialise the helper class for a given particle radius (a) and bend radius (R).
        Both a and R should be read as being relative to half of the duct height,
        i.e. read a as 2a/H and R as 2R/H (H being the duct height).
        The aspect ratio (ar) is with respect to the average duct height
        (or equivalently, the height in the middle) and can either 2 or 4.
        The trapezoid delta/slope (D) is the relative difference in height at
        the right side compared to the centre, it may be 0.0,0.1,0.2,0.3,0.4
        (or also 0.05 when ar=4).
        (Data for additional cross-sections may be added in the future)
        """
        self._a = a
        self._R = R
        self._ar = ar
        self._D = D
        self._load_cs_data()
        self._update_Da_data()
        self._update_R_data()
        self._setup_interpolants()
        return
    def _load_cs_data(self):
        """Load data associate with the desired cross-section"""
        if self._ar==2:
            # If the  data file is not found, a FileNotFoundError exception will be raised
            self._raw_data = np.load('VST_data/strap2x1_lift_data.npz')
        elif self._ar==4:
            # If the data file is not found, a FileNotFoundError exception will be raised
            self._raw_data = np.load('VST_data/strap4x1_lift_data.npz')
        else:
            raise ValueError("The requestied aspect ratio is not available "+\
                             "(only 2 and 4 are currently available)")
        return
    def _update_Da_data(self):
        """Pre-process the raw data associated with the desired D and a"""
        DaR_lookup = self._raw_data['DaR']
        Dis = np.where(DaR_lookup[:,0]==self._D)[0]
        if len(Dis)==0:
            raise ValueError("The requested trapezoid delta/slope is not available "+\
                             "(only 0.0,0.1,0.2,0.3,0.4 are currently available, "+\
                             "or also 0.05 for ar=4)")
        ais = np.where(DaR_lookup[Dis,1]==self._a)[0]
        if len(ais)==0:
            raise ValueError('The requested particle radii is not available for the '+\
                             'specified aspect ratio and trapezoid delta/slope.')
        Ris = Dis[ais] #np.where((DaR_lookup[:,0]==self._D)*(DaR_lookup[:,1]==self._a))[0]
        self._R_values = DaR_lookup[Ris,2]
        Da_data = np.array([self._raw_data['data'][index] for index in Ris])
        # strip NaN data around the edges (there's possibly a cleaner way)
        # This relies on the strap2x1 data being saved with interpolation to fill a rectangular region...
        r = Da_data[0,0,:,20]
        z = Da_data[0,1,20,:]
        rf = np.where(np.isfinite(r))[0]
        zf = np.where(np.isfinite(z))[0]
        self._rs = r[rf[0]:rf[-1]+1]
        self._zs = z[zf[0]:zf[-1]+1]
        self._Da_data = Da_data[:,:,rf[0]:rf[-1]+1,zf[0]:zf[-1]+1]
        # Create a NaN mask
        R = Da_data[0,0,rf[0]:rf[-1]+1,zf[0]:zf[-1]+1]
        Z = Da_data[0,1,rf[0]:rf[-1]+1,zf[0]:zf[-1]+1]
        self._nan_mask = np.ones(Z.shape)
        self._nan_mask[np.abs(Z)>0.975/self._a-1+R*self._D/self._ar] = np.nan
        return
    def _update_R_data(self):
        """Pre-process the raw data associated with the desired R"""
        if self._R<self._ar:
            raise ValueError("The requested bend radius R is non-physical. Choose R>{:d}".format(self.ar_)+
                             "(ideally R>>{:d}) where R is interpreted to be relative".format(self.ar_)+
                             "to half the duct height, i.e. read R as 2R/H")
        if self._R<40.0 or self._R>320.0:
            print("Warning: Extrapolation to bend radii 2R/H<40 or 2R/H>320 may be inaccurate")
        if self._R in self._R_values:
            # use the respective raw data as is
            index = np.where((self._R_values==self._R))[0][0]
            self._DaR_data = self._Da_data[index]
        else:
            # construct an interpolant
            Rs = np.array(self._R_values)
            XS = self._Da_data[0][0]
            ZS = self._Da_data[0][1]
            DaR_data = [XS,ZS]
            for fi in range(2,11):
                F_array = np.array([self._Da_data[Ri][fi] for Ri in range(len(self._Da_data))])
                F_interp = interp1d(1.0/Rs[::-1],F_array[::-1,:,:],axis=0,
                                    kind='cubic' if len(self._R_values)>3 else 'linear',
                                    bounds_error=False,fill_value="extrapolate",assume_sorted=True)
                DaR_data.append(F_interp(1.0/self._R))
            self._DaR_data = np.array(DaR_data)
        return
    def _setup_interpolants(self):
        """Constructs interpolants of the appropriate fields"""
        nhH,nhW = (1.0+self._D-self._a)/self._a,(self._ar-self._a)/self._a
        bbox = [-nhW,nhW,-nhH,nhH]
        self._Cr_RBS = RBS(self._rs,self._zs,self._DaR_data[2,:,:], bbox=bbox)
        self._Cz_RBS = RBS(self._rs,self._zs,self._DaR_data[3,:,:], bbox=bbox)
        self._Lr_RBS = RBS(self._rs,self._zs,self._DaR_data[4,:,:], bbox=bbox)
        self._Lz_RBS = RBS(self._rs,self._zs,self._DaR_data[5,:,:], bbox=bbox)
        self._Sr_RBS = RBS(self._rs,self._zs,self._DaR_data[6,:,:], bbox=bbox)
        self._Sz_RBS = RBS(self._rs,self._zs,self._DaR_data[7,:,:], bbox=bbox)
        self._Up_RBS = RBS(self._rs,self._zs,self._DaR_data[8,:,:], bbox=bbox)
        self._Wr_RBS = RBS(self._rs,self._zs,self._DaR_data[9,:,:], bbox=bbox)
        self._Wz_RBS = RBS(self._rs,self._zs,self._DaR_data[10,:,:],bbox=bbox)
        self._kappa = 4.0/(self._R*self._a**3)
        self._Fr_RBS = RBS(self._rs,self._zs,self._DaR_data[4,:,:]+self._kappa*self._DaR_data[6,:,:],bbox=bbox)
        self._Fz_RBS = RBS(self._rs,self._zs,self._DaR_data[5,:,:]+self._kappa*self._DaR_data[7,:,:],bbox=bbox)
        return
    def _plot_DaR_component(self,ind,apply_mask=True):
        """Useful for testing/debugging purposes"""
        RS,ZS = np.meshgrid(self._rs,self._zs,indexing='ij')
        plt.figure(figsize=(1+3*self._ar,3))
        if apply_mask:
            plt.contourf(RS,ZS,self._DaR_data[ind,:,:]*self._nan_mask,17)
        else:
            plt.contourf(RS,ZS,self._DaR_data[ind,:,:],17)
        plt.colorbar()
        nhHL,nhHR,nhW = (1.0-self._D)/self._a,(1.0+self._D)/self._a,self._ar/self._a
        plt.plot([-nhW+1,nhW-1,nhW-1,-nhW+1,-nhW+1],[-nhHL+1,-nhHR+1,nhHR-1,nhHL-1,-nhHL+1],'r--')
        plt.plot([-nhW,nhW,nhW,-nhW,-nhW],[-nhHL,-nhHR,nhHR,nhHL,-nhHL],'k-',lw=1)
        plt.xlim(-nhW,nhW)
        plt.ylim(-nhHR,nhHR)
        plt.gca().set_aspect(1.0)
        plt.show()
    def get_cross_section(self):
        """Get the current cross-section details
        (returns the aspect ratio and trapezoid delta/slope parameters)"""
        return self._ar,self._D
    def get_bounds(self): # TODO: probably want to re-think this...
        """Get the bounds of the current cross-section
        (given in the form [r_min,r_max,z_min,z_max])"""
        return [-self._ar/self._a,self._ar/self._a,-1.0/self._a,1.0/self._a]
    def get_particle_bounds(self): # TODO: probably want to re-think this...
        """Get the bounds for the particle centre in the current
        cross-section (given in the form [r_min,r_max,z_min,z_max])"""
        return [-self._ar/self._a+1.0,self._ar/self._a-1.0,-1.0/self._a+1.0,1.0/self._a-1.0]
    def check_particle_bounds(self,x,z):
        """Check if the given particle coordinates (of its centre)
        are physically reasonable for the current cross-section."""
        xmax = self._ar/self._a-1
        S = self._D/self._ar
        zmax = 1.0/self._a-1.0+x*S
        #return (np.abs(x)<xmax) and (np.abs(z)<=zmax) # a bit crude (i.e. not entirely correct)
        return (np.abs(x)<xmax) and (zmax-np.abs(z)>=(1.0+S**2)**0.5) # better
    def get_aspect_ratio(self):
        """Get the current (mean) aspect ratio (ar)."""
        return self._ar
    def set_aspect_ratio(self,ar,D=None,a=None,R=None):
        """Change the aspect ratio.
        Optionally the D,a,R parameters may be simultaneously updated."""
        self._ar = ar
        self._load_cs_data()
        if D is not None:
            self._D = D
        if a is not None:
            self._a = a
        self._update_Da_data()
        if R is not None:
            self._R = R
        self._update_R_data()
        self._setup_interpolants()
        return
    def get_D(self):
        """Get the current trapezoid delta/slope parameter (D)."""
        return self._D
    def set_D(self,D,a=None,R=None):
        """Change the trapezoid delta/slope parameter (D).
        Optionally the a,R parametters may be simultaneously updated."""
        self._D = D
        if a is not None:
            self._a = a
        self._update_Da_data()
        if R is not None:
            self._R = R
        self._update_R_data()
        self._setup_interpolants()
        return
    def get_available_a(self): # TODO: possibly rethink this...
        """Get a list of available particle radii"""
        return np.unique(self._raw_data['DaR'][:,1])
    def get_a(self):
        """Get the current particle radius"""
        return self._a
    def set_a(self,a,R=None):
        """Change the particle radius.
        Optionally update R simultaneously"""
        self._a = a
        self._update_Da_data()
        if R is not None:
            self._R = R
        self._update_R_data()
        self._setup_interpolants()
        return
    def get_R(self):
        """Get the current bend radius"""
        return self._R
    def set_R(self,R):
        """Change the bend radius of the duct (and update the interpolants)"""
        self._R = R
        self._update_R_data()
        self._setup_interpolants()
    def migration_force(self,r,z):
        """Get the (net) migration force for a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via rho U_m^2 a^4 / H^2 )"""
        return np.squeeze([self._Fr_RBS(r,z),self._Fz_RBS(r,z)])
    def migration_force_jacobian(self,r,z):
        """Get the jacobian of the (net) migration force for a neutrally
        buoyant spherical particle centred at (r,z) within the
        cross-section (non-dimensionalised via rho U_m^2 a^3 / H^2 )"""
        return np.squeeze([[self._Fr_RBS(r,z,dx=1),self._Fr_RBS(r,z,dy=1)],
                           [self._Fz_RBS(r,z,dx=1),self._Fz_RBS(r,z,dy=1)]])
    def migration_velocity(self,r,z):
        """Get the migration velocity for a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m a / H )"""
        return np.squeeze([self._Fr_RBS(r,z)/self._Cr_RBS(r,z),
                           self._Fz_RBS(r,z)/self._Cz_RBS(r,z)])
    def drag_coefficient(self,r,z):
        """Get the drag coefficients in the r,z directions of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via mu a )"""
        return np.squeeze([self._Cr_RBS(r,z),self._Cz_RBS(r,z)])
    def secondary_flow_drag(self,r,z):
        """Get the drag coefficients in the r,z directions of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via rho U_m^2 a^4 / H^2 )"""
        return np.squeeze([self._Sr_RBS(r,z),self._Sz_RBS(r,z)])
    def axial_velocity(self,r,z):
        """Get the terminal/steady axial velocity of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m a / H )"""
        return np.squeeze(self._Up_RBS(r,z))
    def spin_components(self,r,z):
        """Get the terminal/steady r,z spin components of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m / H )"""
        return np.squeeze([self._Wr_RBS(r,z),self._Wz_RBS(r,z)])
    def plot_migration_force(self,apply_mask=True,show=True):
        """Produces a rough sketch of the magnitude of the migration force including
        the zero contours of the horizontal and vertical components."""
        Fr = self._DaR_data[4,:,:]+self._kappa*self._DaR_data[6,:,:]
        Fz = self._DaR_data[5,:,:]+self._kappa*self._DaR_data[7,:,:]
        if apply_mask:
            Fr *= self._nan_mask
            Fz *= self._nan_mask
        RS,ZS = np.meshgrid(self._rs,self._zs,indexing='ij')
        plt.figure(figsize=(1+3*self._ar,3))
        plt.contourf(RS,ZS,(Fr**2+Fz**2)**0.5,17)
        plt.colorbar()
        plt.contour(RS,ZS,Fr,[0.0],colors=['k'])
        plt.contour(RS,ZS,Fz,[0.0],colors=['w'])
        nhHL,nhHR,nhW = (1.0-self._D)/self._a,(1.0+self._D)/self._a,self._ar/self._a
        plt.plot([-nhW+1,nhW-1,nhW-1,-nhW+1,-nhW+1],[-nhHL+1,-nhHR+1,nhHR-1,nhHL-1,-nhHL+1],'r--')
        plt.plot([-nhW,nhW,nhW,-nhW,-nhW],[-nhHL,-nhHR,nhHR,nhHL,-nhHL],'k-',lw=1)
        plt.xlim(-nhW,nhW)
        plt.ylim(-nhHR,nhHR)
        plt.gca().set_aspect(1.0)
        if show:
            plt.show()
        return
    # end of class
