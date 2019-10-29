#function to compute antenna response pattern and voltages

import numpy as np
from scipy.integrate import quad as qd
import scipy.special as sp


def response_pattern(theta,D,foclen,freq,dy,dx, Z_0 = 377,c= 3e8 ):
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

    power =[]
    m = theta
    power_pattern =[]
    for m_i in range(m.size):
         tot_electric_int = qd(real_int,-D/2.0,D/2.0,args=(m[m_i],dy,dx))[0] + 1j*qd(imag_int,-D/2.0,D/2.0,args=(m[m_i],dy,dx))[0]
        
         norm =  qd(real_int,-D/2.0,D/2.0,args=(0.0,0.0,0.0))[0] + 1j*qd(imag_int,-D/2.0,D/2.0,args=(0.0,0.0,0.0))[0]
         power_pattern.append((np.conj(tot_electric_int)*tot_electric_int)/(np.abs(norm)**2)*Z_0)
        
    return [tot_electric_int, nom_power_pattern]


def get_effe_len(E_field,Z_ants,Z_0 = 377,k= 1.38064852e-23):
	"""
	This function compute the effective lenght of single polaraztion half-wavelength
        
        Parameters:
	-----------
	E_field: array
		Electric as function zenith angle at observing frequency f.

	Z_ant : integer
             Antenna impendancein unit of Ohms.

	Z_0 : integer
             Impendance of free space in unit of Ohms.

   	k: float
	      Boltzmann constant in unit of meter square times kg per second square per Kelvin.
       

	
	returns:
	eff_len: array
		The effective lenght as function of zenith angle at the observing frequency f.


	"""

	#Current on dipole conductors induced by the Electric at far field E_field.
	Z_0 = 377 #Ohms
      
	E_amp_in = np.abs(E_field)

	I_p = 2.0*E_amp_in*np.sqrt(Z_0)/(Z_ant +Z_0)
	
	#effective length
	eff_len = -(E_far*4.0*np.pi)/(k*Z_0*I_p)
	
	return eff_len


def get_ants_voltage(E_field,eff_len,z_angle_arr):
	"""
	This function compute the antenna voltage by taking the dot product of the Electric Field 		(as function of zenith angle) and half-wavelength dipole effective length (as function of 		zenith angle)

        
        Parameters:
	-----------
	E_field: array
		Electric as function zenith angle at observing frequency f.

	eff_len : array
             Effective lenght as fuction of zenith angle at observing frequency f.
	z_angle_arr: array
	     The zenith angle array in radians from -pi/2 to pi/2.
	
	returns:
	ant_voltage: array
		The antenna voltage output of single polarazation as function of zenith angle at the 		observing frequency f.
	ant_tot_voltage_out: complex 
		The total antenna voltage output of single polarazation as function of zenith 			angle at the observing frequency f.

	


	"""
	if E_field.shape == eff_len.shape:
		pass
	else:
		print "E_field and eff_len have different shape"
		
		break


	ant_voltage_outp = np.dot(E_field,eff_len)
	ant_tot_voltage_outp= sum( 2.0*np.pi*ant_voltage_outp*np.sin(z_angle))

	return [ant_voltage_outp,ant_tot_voltage_outp]

def get_visibility(ant_tot_voltage_outp_i,ant_voltage_j,bl_ij,freq,z_angle_arr):
	"""
	This function compute the visiblity from the pair of participating antenna's voltages by 		taking the complex product of the total voltage out of antenna i with antenna j, multiply the 		product  phase factor due to geometrical delay of the electric field coming from 		cosmic source and with antenna closely packed on the EW NS plane. assuming symetric about
	zenith vector, multiply by the solid angle 2pisin(theta) and summing over the -pi/2, pi/2. 		
        
        Parameters:
	-----------
	ant_tot_voltage_outp_i: complex
		The total antenna voltage output of single polarazation as function of zenith 			angle at the observing frequency f.
	ant_tot_voltage_outp_j: complex
		The total antenna voltage output of single polarazation as function of zenith 			angle at the observing frequency f.

	bl_ij : array shape(1,3)
          	The distance between antenna_i position and antenna_j in topocentric cordiante 			system, meters.
	freq: float
		Observing frequency in MHz.
	z_angle_arr: array
	     The zenith angle array in radians from -pi/2 to pi/2.


	
	returns:
	ant_voltage: array
		The antenna voltage of single polarazation as function of zenith angle at the 		observing frequency f.

	

	"""
	
	vis = np.sum(2.0*np.pi*np.conj(ant_tot_voltage_outp_i)*ant_tot_voltage_outp_j*np.exp(1j*2.0*np.pi*np.abs(bl_ij)*np.cos(z_angle_arr))*np.sin(z_angle))
	
	return vis








	


















