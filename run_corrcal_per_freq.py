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
    	Cov_dic : array; shape(NFREQS,NANTS)
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
    
	
def remove_phase_degen(gain_sols,ants_xyz):
	"""
	This fuction remove the phase gradient.
	Method   :  fits the predicted phase from the Array layout
	

     	Parameters
    	----------
    	gain_sols      : array; shape:(NFREQ,NANTS)
		       Antenna gain solutions for NFREQS frequency channels.

	ants_xyz       : array; shape:(NANTS,3)
	      Antenna positions in topocentric cordinates (meters).

	
    
	Returns
    	-------
    	 phase_sols : array; shape(NFREQS,NANTS)
        	Antenna gain phase solution.

  
	"""



	nant= gain_sols.shape[1]
	nfreq = gain_sols.shape[0]
	phase=np.imag(gain_sols)

	lhs=np.dot(ants_xyz.transpose(),ants_xyz)

	pp=0*phase
	phase_sols=[]
	for freq_i in range(nfreq):
    		rhs=np.dot(ant_xyz.transpose(),phase[freq_i,:])
    		fitp=np.dot(np.linalg.pinv(lhs),rhs)
    		pred=np.dot(ant_xyz,fitp) 
    		pp[freq_i,:]= phase[freq_i,:]-pred #removing phase gradient
                phase_sol.append(pp[freq_i,:])


	return phase_sols





def get_gain_autocorr_per_ant(gain_sols):


	Gain={}
	ant = np.arange(gain_sols.shape[1])
	for ii in range(len(ant)):
        	tmp=[]
        	for freq_i in range(gain_sols.shape[0]):
            		for i in range(gain_sols.shape[1])):
                		if ii == i:
                    		tmp.append(gain_sols[freq_i][i])
        	Gain[ii]= tmp


	GAIN_SOL ={}

	for ant_i in range(len(ant)):
		tmp = []
		for key in Gain.keys():
			if ant_i == key:
				tmp.append(Gain[key])
		

		GAIN_SOL[ant_i]=tmp

	Autocorr_amp={}
	Autocorr_phase{}
  
 	for ant_j in GAIN_SOL.keys():
		
        	for k in range(GAIN_SOL[ant_j].shape[0]-1):
			
			phase = np.angel(GAIN_SOL[ant_j])
			phase = phase - np.mean(phase)
			amp = np.absolute(GAIN_SOL[ant_j])
			amp = amp - np.mean(amp)
			tmp0=[]
			tmp1 =[]

			for n in range(n_data-k):
				tmp0.append(phase[n]*phase[n+k])
				tmp1.append(amp[n]*amp[n+k])
		
				
			phase_autocorr[k]=np.sum(tmp0)/(n_data-k)
			amp_autocorr[k]=np.sum(tmp1)/(n_data-k)

		Autocorr_amp[ant_j]   =  amp_autocorr
		Autocorr_phase[ant_j] =  phase_autocorr

	

	return Autocorr_amp,Autocorr_phase

	

def get_gain_autocorr_all_ants(gain_sols):








	
for ii in range(len(ant)):
	tmp = []
	for key in Gain.keys():
			if ii == key:
				tmp.append(Gain[key])
	ant_ii = tmp # np.column_stack(tmp)
	#print ant_ii.shape, ii

	GAIN_SOL[ii]=ant_ii



N_autocorr_amp={}
N_autocorr_phase={}
for ant_i in range(len(GAIN_SOL)):
	Autocorr_amp=[]
	Autocorr_phase =[]
	
	for ni in range(len(GAIN_SOL[ant_i])):
		GAIN_SOL_AVG = GAIN_SOL[ant_i][ni]/np.mean(GAIN_SOL[ant_i][ni])
        #print GAIN_SOL_AVG
        phase_autocorr ={}
        amp_autocorr ={}
        n_data = len(GAIN_SOL_AVG)
        for k in range(n_data-1):
			#phase = np.arctan2(np.real(GAIN_SOL_AVG[ant]),np.imag(GAIN_SOL_AVG[ant]))
			phase = np.real(GAIN_SOL_AVG) 
			phase = phase - np.mean(phase)
			#amp = np.sqrt(np.real(GAIN_SOL_AVG[ant])**2 + np.imag(GAIN_SOL_AVG[ant])**2)
			amp = np.real(GAIN_SOL_AVG)
			amp = amp - np.mean(amp)
			tmp0=[]
			tmp1 =[]
			for n in range(n_data-k):
				tmp0.append(phase[n]*phase[n+k])
				tmp1.append(amp[n]*amp[n+k])
		
				
			phase_autocorr[k]=np.sum(tmp0)/(n_data-k)
			amp_autocorr[k]=np.sum(tmp1)/(n_data-k)
        Autocorr_amp.append(amp_autocorr)
        Autocorr_phase.append(phase_autocorr)
	N_autocorr_amp[ant_i]= Autocorr_amp
	N_autocorr_phase[ant_i]= Autocorr_phase	
