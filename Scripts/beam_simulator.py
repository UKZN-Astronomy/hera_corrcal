import numpy as np
from scipy.integrate import quad as qd
from scipy import special as sp


#2dGuassian beam

#primary beam function

def gaussian_beam(theta_c, phi_c,theta,phi,sigma_fwhm,e_theta,e_phi,e_sigma):
    """
    This function compute the 2d gassian beam for full width half maximum,sigma_fwhm,
    value at location theta,phi (radians)
    """
    
    #Guassian Primary beam
    gpbeam = np.exp(-((theta -theta_c + e_theta)**2 + (phi - phi_c+ e_phi)**2)/(2.0*(sigma_fwhm + e_sigma)**2))
    
    return gpbeam

def sigma_func(d,freq_i):
    """
    This function compute the sigma_fwhm base on the Limits of Resolution: The Rayleigh Criterion
    
    """
    lambda_ = (3.0*10**8)/freq_i # frequency in Hertz
    sigma = 0.44*(lambda_/d)
    return sigma
    


#parametized airydisk function

def param_airydisk(theta_c,phi_c,theta_x,phi_y,a_x,a_y,d_diam,d_block,freq_i, theta_error =0.0,phi_error=0.0):
	"""
	This function compute the airy disk function given the angular coordinate theta_x and phi_y in radians, center coodinate 		(theta_c,phi_c),  major and semi-major parameters a_x and a_y in meters, dish diameter d_diam, diameterof blocking cage d_block  		frequency freq_i
	

	"""

    
    	wavelen = 3e8/freq_i
    	b= d_block/d_diam
    	k = (2.0*np.pi)/wavelen
    	sqrt_term = np.sqrt((a_x**2)*(theta_x- (theta_c + theta_error) )**2 + (a_y**2)*(phi_y - (phi_c +phi_error))**2)
    
    	power_pattern = np.power(2.0*(sp.j1(k*sqrt_term)-b*sp.j1(b*k*sqrt_term))/((1-b**2)*k*sqrt_term),2)
    
    	return power_pattern



#geometrical optics beam model
def complex_modulus(z):
    
    a= np.real(z)
    b= np.imag(z)
    
    return np.sqrt(a**2 + b**2)




def geomtrical_optics_beam_model(theta,D,foclen,freq,dy,dx, Z_0 = 377,c= 3e8 ):
    """
    This fuction compute the dish response pattern given the co-altitude angle 
    theta (zenith angele, dish diameter D, observing frequency in MHz, and the axial dy
    and letarial dx feed displacement
    

     Parameters
    ----------
    theta : 1d array,
        Arrays of zenith angle in radians. 
    D: float
	Dish diameter in meters.
    
    foclen: float
           Dish focus lenght in meters.

    freq: float
	 
	Observing frequency in MHz.

   
   dy: float
       The feed axial dispacement in y-direction.

   dx:float
           The feed lateral dispacement in x-direction.

   Z_0 : integer
       Impendance of free space in unit of Ohms
   c : float
	Speed of light in vacuum in unit of meters/seconds
       

         
         
    
    
    Returns
    -------
    tot_electric_int : array
        Array with total Electric as function of zenith angle.

    nom_power_pattern: array
	 Array with normalized antenna power pattern as function of zenith angle.

    """

   
    lambda_= c/(freq*10**6)  #meters

    #farfield distance

    d_f = 2.0*D**2/lambda_
    y_0 = foclen
    y = lambda x : x**2/4.0/y_0 #dish surface
    
    #focus dispacement
    df = lambda dy,dx :np.sqrt((y_0 +dy)**2 + dx**2)
    #total path travelled by a planewave from a source at far-field to the focal point
    tot_path = lambda x,theta,dy,dx : abs(np.sin(theta)*x - np.cos(theta)*y(x) + h) + np.sqrt(x**2 + (df(dy,dx)- y(x))**2)

    exp_ = lambda x,theta,dy,dx : np.exp(1j*2.0*(np.pi/lambda_)*tot_path(x,theta,dy,dx))

    real_int = lambda x,theta,dy,dx : exp_(x,theta,dy,dx).real*np.sqrt(1.0 + (x/2./y_0)**2)
    imag_int = lambda x,theta,dy,dx : exp_(x,theta,dy,dx).imag*np.sqrt(1.0 + (x/2./y_0)**2)

    h= d_f
 
    tot_electric_int = qd(real_int,-D/2.0,D/2.0,args=(theta,dy,dx))[0] + 1j*qd(imag_int,-D/2.0,D/2.0,args=(theta,dy,dx))[0]
    power_pattern = complex_modulus(tot_electric_int)**2/Z_0
        
    norm =  complex_modulus(qd(real_int,-D/2.0,D/2.0,args=(0.0,0.0,0.0))[0] + 1j*qd(imag_int,-D/2.0,D/2.0,args=(0.0,0.0,0.0))[0])**2/Z_0
    norm_power_pattern= power_pattern/norm
        
    return norm_power_pattern




                      

    
    
