#This scipt compute the expected sky covariance in spherical harmonics space
#computing expected sky covaraince matrix
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import healpy as hp
import time





def get_expected_sky_cov(ublDict,beammodel,antspos,long_0,lalt_0,cl_21cm,cl_gal_syntron,cl_gal_ff,cl_extragal_ps,cl_extragal_ff,freq,nside):
    """
    This fuction compute the sky expected covariance given the redundant block information,
    antenna beam models per frequency antenna position, array location (longitude,latitude), 
    the expect angular power specttrum  during the time observations:cl_21cm,
    cl_gal_syntron,cl_gal_ff, cl_extragal_ps,cl_extragal_ff at frequency    

     Parameters
    ----------
    ublDict: dictionary (NBlocks,Nbaseline,baselinevec:(ant1_index,ant2_index)
       
    beammodel : array (NANTS,NBEAMSIZE)
	      Beam model for each NANTS, with beam array of NBEAMSIZE sky position, this is zenith 
	      angel.
    antspos: array (NANTS,3)
	      Antenna positions in topocentric cordinate (meters).
    
    long_0,lalt_0 : floats
	      Array location, (longitude,latitude).

    cl_21cm: array
	    Expected 21 cm angular power spectrum.

    cl_gal_syntron: array
	    Expected galactic syntron emission angular power spectrum.

    cl_gal_ff: array
	    Expected galactic free-free emission angular power spectrum.

    cl_extragal_ps: array
	    Expected extra-galactic radio point sources angular power spectrum.

    cl_extragal_ff: array
	    Expected extra-galactic free-free emission angular power spectrum.
  
    freq : float
	 Observing frequency

    nside: int
	 nside for resolution in computing spherical hamonic coefficeint of interferometer
	 beams.   
    

    
    Returns
    -------
    Cov_dic : Dictionary {redbaseline:array(NBLOCK,NBLOCK)}
        A dictionary of covariance matrix per redundant baseline.

    I_k     : Dictionary {redbaseline:array(NBLOCK,NBLOCK)}
	  A dictionary of correlation matrix per redundant baseline.

    """


    
    
    t1 =time.time() 
    #observing wavelenght
    
    lambda_ = (3e8)/(freq*10**6)
    
    #creating sky position for given nside, which should match nside of beammap healpix
    lmax = 3*nside**2 -1
    npix = 12*nside**2
    s= np.array(hp.pix2vec(nside,np.arange(npix)))
    
    
    Cov_dic ={} # dictionary of expected of the sky covariance matrices for all redundant blocks at observing frequency
    I_dic ={} # dictionary of expected of the correlation ...
    
    
    for ubl_k, block_k in ublDict.iteritems(): 
        cov_k =np.zeros((len(block_k),len(block_k)), dtype='complex') #covarinace matrix for a particular redundant block
        I_k =np.zeros((len(block_k),len(block_k)), dtype='complex') # correlation matrix ...
        for bl_w in range(len(block_k)):
             for bl_z in range(len(block_k)):
                    
                    
                    #basiline alpha
                    delta_alpha = antspos[block_k[bl_w][0]] - antspos[block_k[bl_w][1]] # baseline for interferemeter pair alpha
                    bl_alpha = np.resize(np.repeat(delta_alpha,npix),(3,npix))
                    fringe_alpha =np.exp(2.0*np.pi*1j*np.sum(bl_alpha*s)/lambda_) # fringe for bl_alpha
                   
                    beammaprot_i = get_rotmap(beammodel[block_k[bl_w][0]],[long_0,lalt_0]) #rotated beam for antenna i to array location on earth longitude and latitude
                    beammaprot_j = get_rotmap(beammodel[block_k[bl_w][1]],[long_0,lalt_0]) #rotated beam for antenna j to array location on earth longitude and latitude
                    interfer_beam_alpha = np.conj(beammaprot_i)*beammaprot_j*fringe_alpha # interferometry beam for pair of antenna (i,j)
                    
                    blm_alpha = hp.map2alm(interfer_beam_alpha,lmax=lmax,mmax=0,iter=3) #spherical harmonic transform of interferomtric beam

                    #baseline beta
                    delta_beta= antspos[block_k[bl_z][0]] - antspos[block_k[bl_z][1]] # baseline for interferemeter pair beta
                    bl_beta = np.resize(np.repeat(delta_beta,npix),(3,npix))
                    fringe_beta =np.exp(2.0*np.pi*1j*np.sum(bl_beta*s)/lambda_)
                    
                    beammaprot_n = get_rotmap(beammodel[block_k[bl_z][0]],[long_0,lalt_0])
                    beammaprot_m = get_rotmap(beammodel[block_k[bl_z][1]],[long_0,lalt_0])
                    interfer_beam_beta = np.conj(beammaprot_n)*beammaprot_m*fringe_beta
                    
                    blm_beta = hp.map2alm(interfer_beam_beta,lmax=lmax,mmax=0,iter=3)

                   
                    #computing covariance and correlation in spherical harmonics space
                    cl_sky =  cl_extragal_ps #[cl_extragal_ff[l_j] + cl_extragal_ps[l_j] for l_j in range(len(cl_extragal_ff))] # sky power spectrum with galactic centre below the horizon and ignoring 21 cm signal
		    #print 'cl', len(cl_sky), len(blm_alpha[1:]), len(np.conj(blm_alpha[1:])*blm_beta[1:])
                    tmp_alpha_beta = np.conj(blm_alpha[1:])*blm_beta[1:]*cl_sky #computing the covariance of the sky in SH-Space
                    cov_alpha_beta = np.sum(tmp_alpha_beta)
                    tmp_alpha_alpha = np.conj(blm_alpha[1:])*blm_alpha[1:]*cl_sky
                    cov_alpha_alpha= np.sum(tmp_alpha_alpha)
                    tmp_beta_beta = np.conj(blm_beta[1:])*blm_beta[1:]*cl_sky
                    cov_beta_beta = np.sum(tmp_beta_beta)
                    I_alpha_beta = cov_alpha_beta/np.sqrt(cov_alpha_alpha*cov_beta_beta)
                    print  'cov_alpha_beta', cov_alpha_beta,'cov_alpha_alpha',cov_alpha_alpha,'cov_beta_beta',cov_beta_beta, ubl_k
                    cov_k[bl_w][bl_z] = cov_alpha_alpha
                    I_k[bl_z][bl_w] = I_alpha_beta
                
                    
        Cov_dic[ubl_k]= cov_k
        I_dic[ubl_k] = I_k
    t2 =time.time()  
    print 'time elapse ',  t2-t1
    return Cov_dic,I_dic











def get_redblock_inf(redblock_cov_matrix_eig,threshold):
    """
    This fuction reduce the size of redundant block covariance matrix by select the only eigenmodes
    with magnitude above expected signal threshold and set to zero the rest of modes.
    

     Parameters
    ----------
    redblock_cov_matrix_eig: array(NFREQS,NBLOCKS,NEIGVALS,NEIGVECS) dtype: complex
        Arrays of NFREQS, NBLOCK of redundant blocks, NEIGVALS & NEIGVECS for each redundant block.
    threshold              : float
			    Threshold of the eigenmode that will be selected.
   
         
         
    
    
    Returns
    -------
    R_all_freq : (NFREQS,NBLOCKS,2*NEIGVALS,2*NEIGVECS) array
        Array for all frequency with redundant  block reduced in size by selecting nmax of largest
	eigenmodes. The real and imaginary has be separated, hence 2 times the size of eigmodes.
    edges: array 
	 Array containing the index of the block edges.

    """


    
    RR= redblock_cov_matrix_eig
    nblock=RR[0].shape[0]

    edges=np.zeros(nblock+1,dtype='int')
    R_all_freq=[]
    for nu_i in range(len(RR)):
    	nmax=0
    	thresh=threshold
    	nvis=0
    	for i in range(nblock):
        	myeig=qq[nu_i][i][0]
        	nkept=np.sum(myeig>thresh*myeig.max())
        	nvis=nvis+myeig.size
        	edges[i+1]=edges[i]+myeig.size
        	if nkept>nmax:
            		nmax=nkept
    	print 'total number of visibilities is ' + repr(nvis) + ' with ' + repr(nmax) + ' kept eigenvalues at max.'
    	# computeting R = vl^1/2, l_above_threshold
    	R=np.zeros([2*nvis,nmax*2])
    
   	myeig_vec =[]
    	for i in range(nblock):
        	myeig=np.real(RR[nu_i][i][0])
        	myvecs=np.real(RR[nu_i][i][1])
        	ind=myeig>thresh*myeig.max() # picking up an index
        	#print 'index max eiegen', ind
        	myeig_use=myeig[ind] # pict max eigenvalue
        	#print myeig_use.shape
        	myvecs_use=myvecs[:,ind] # pick  corresponding max eigenvec
        	#print myvecs_use.shape
        	for j in range(len(myeig_use)):
            		myvecs_use[:,j]= myvecs_use[:,j]*np.sqrt(myeig_use[j])
            
            		q[2*edges[i]:2*edges[i+1]:2,2*j]=  np.column_stack(myvecs_use[:,j])
            		q[(2*edges[i]+1):2*edges[i+1]:2,2*j+1]= np.column_stack(myvecs_use[:,j])

	R_all_freq.append(R)


    return R_all_freq, edges




    

