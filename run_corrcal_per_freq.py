#This function compute antenna gain solution per frequency using correlation calibration (corrcal)

from scipy.optimize import fmin_cg
import numpy as np
import matplotlib.pyplot as plt
import corrcal2 as corrcal
from scipy import optimize
import time


def get_corrcal_gainsol(data,ant1,ant2,gain_int,noise,sky_cov_vecs,src_vecs,block_edges,gain_fac,maxiter=1000):
	"""
    	This fuction solve for antenna gain per frequency.
	Method   :  Conjugate-gradient Method 
	Software :  python scipy.optimize --- fmin_cg

     	Parameters
    	----------
    	data      : array; shape:(NFREQ,2*NVIS)
		   Visibility data for NFREQ frequency channnels, and grouped into redundant blocks and seperated in real and imaginary part.

	ant1 : array; shape:(NVIS,1)
	      Antenna indeces for ant1 set sorted according to redundant blocks.

	ant2 : array; shape:(NVIS,1)
	      Antenna indeces for ant2 set sorted according to redundant blocks
       
    	gain_int  : array; shape(NFREQS,2*NANTS)
	      Antenna gain inintial guess for NFREQS frequency channels, and with real and imaginary part separated.

    	noise    :  array; shape:(NFREQ,2*NVIS)
			Per Visibility noise for NFREQ frequency channnels, and grouped into redundant blocks and seperated in real and imaginary part.
			
	     
    
	sky_cov_vecs : array; shape:(NFREQ,2*NVIS)

	      Expected sky data  from array covariance matrix for NFREQ frequency channnels, and grouped into redundant blocks and seperated in real and imaginary part.
	src_vecs : array; shape:(NFREQ,2*NVIS)

	      Expect sky data  from bright point sources for NFREQ frequency channnels, and grouped into redundant blocks and seperated in real and imaginary part.

	edges : array; shape:(NVIS,1)
	      The edges of the redundant blocks.

	gain_fac :int
		Scaling factor for gain so that they take appropriate steps.

	maxiter : int
		Maximum iteration for the steps that fmin_cg should before the optimum gain solution.
		

	     


    
	Returns
    	-------
    	Cov_dic : array; shape(NFREQS,2*NANTS)
        	Antenna gain solution

  
	"""

	nants = gain_int.shape[1]/2
	gain_sol_cg =[]
	nfreq = data.shape[0]
	for freq_i in range(nfreq):


		#computing the effective noise covariance matrix
 		noise_eff =corrcal.sparse_2level(noise[freq_i),sky_covvecs[freq_i].transpose().copy(),src_vecs[freq_i],2*block_edges)
		
		fac=fac   #scaling factor for gain so that they take appropriate steps
    		normfac=1
    		t1=time.time()
    		gain_sol=fmin_cg(corrcal.get_chisq,gvec*fac,corrcal.get_gradient,(data,noise_eff,ant1,ant2,fac,normfac),maxiter=maxiter)
    		t2=time.time()
   
    		print 'elapsed time to do nonlinear fit for ' + repr(nant) + ' antennas was ' + repr(t2-t1)
    		gain_sol=gain_sol/fac
    		gain_sols=gain_sol[0::2]+np.complex(0,1)*gain_sol[1::2]

		#removing the degeneracies by dividing by an average of gains
									
		gain_sols= gain_sols/np.absolute(np.average(gain_sols))  
		#gain solution


		gain_sol_cg.append(gain_sols)
		

	return gain_sol_cg
    
	

