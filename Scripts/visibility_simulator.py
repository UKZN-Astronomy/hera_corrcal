import numpy as np
import beam_simulator as bs




#Visibility simulation with Gaussain beam

def vis_sim_gaussain_beam(ublDict,Ant_pos,dish_diam,RA_0,DEC_0,RA,DEC,flux,freq,theta_error,phi_error,beam_size_error):
    
        sky_pos_r = [[np.cos(lalt_i)*np.cos(long_i), np.cos(lalt_i)*np.sin(long_i), np.sin(lalt_i)] for (lalt_i,long_i) in zip(DEC,RA)]
        Ant_index = np.arange(len(Ant_pos))
        #transforming antenna position from ENU to equatorial system XYZ
        lat = DEC_0
        XYZ = np.array([XYZ_trans_mat(Ant_pos[i],lat)[0] for i in range(len(Ant_pos))])
        
        #looping over redundant blocks
        vis_array ={}
        for redblock, ants_redblock in ublDict.iteritems():
            
            #computing visibility per baseline within redundant block
            vis_redblock=[]
            for bl_i in range(len(ants_redblock)):
                
                #computing the baseline and projecting it to the phase center\
                ant_1 = ants_redblock[bl_i][0]
                ant_2 = ants_redblock[bl_i][1]
                bl_xyz = XYZ[ant_2]- XYZ[ant_1]
                bl_xyz_proj = baseline_proj(bl_xyz,RA_0,DEC_0)
                
                
                #compute the visibility for all frequency
                vis =np.zeros(len(freq),dtype='complex')
                for freq_i in range(len(freq)):
                    
                    lambda_ = (3e8)/(freq[freq_i]*10**6)
                    sigma_fwhm = bs.sigma_func(dish_diam,freq[freq_i])
                    theta, phi = DEC, RA
                    
                    phase = np.array([np.exp(1j*2.0*np.pi*np.dot(bl_xyz_proj[0],sky_pos_r[r_i]/np.linalg.norm(sky_pos_r[r_i]))/lambda_)  for r_i in range(len(sky_pos_r))])
                    perceive_sky = np.array([flux[i][freq_i]*bs.gaussian_beam(theta[i],phi[i],sigma_fwhm,theta_error[ant_1],phi_error[ant_1],beam_size_error[ant_1])*bs.gaussian_beam(theta[i],phi[i],sigma_fwhm,theta_error[ant_2],phi_error[ant_2],beam_size_error[ant_2]) for i in range(theta.size)])
                    vis[freq_i] = np.sum(perceive_sky*phase)
                    
                vis_redblock.append(vis)
            
            vis_array[redblock] = vis_redblock
            
            
            
        return vis_array
    
HA = lambda lst,ra_obj: lst - ra_obj

def vis_sim_gaussian_beam_bl(ant1_pos,ant2_pos,dish_diam,RA_0,DEC_0,RA,DEC,src_flux_spec_index,freqs,lsts,ant1_theta_error,ant1_phi_error, ant2_theta_error,ant2_phi_error,Ant1_beam_size_error,Ant2_beam_size_error):
    
        vis =np.zeros((lsts.size,freqs.size),dtype='complex')
        
        for lst_i in range(lsts.size):
            
             #star hour angle
            ra=[lsts[lst_i] -RA[src_i] for src_i in range(RA.size)]
            
            src_pos_lmn = [pos_lmn(ra[src_i],DEC[src_i],np.deg2rad(RA_0),np.deg2rad(DEC_0)) for src_i in range(len(ra))]
            
            
            #transforming antenna position from ENU to equatorial system XYZ
            lat = np.deg2rad(DEC_0)
            xyz_1,xyz_2 = XYZ_trans_mat(ant1_pos,lat)[0],XYZ_trans_mat(ant2_pos,lat)[0]
            
            #computing the baseline
            bl_xyz = xyz_2 - xyz_1

           
            #compute the visibility for all frequency
            vis_freq = np.zeros(freqs.size,dtype ='complex')
            for freq_i in range(len(freqs)):
		lambda_ = (3e8)/(freqs[freq_i]*10**6)
                sigma_fwhm = bs.sigma_func(dish_diam,freqs[freq_i]*10**6)
                
                vis_freq[freq_i]= np.sum([flux_func(freqs[freq_i],src_flux_spec_index[src_i][0], src_flux_spec_index[src_i][1])*bs.gaussian_beam(np.deg2rad(DEC_0),np.deg2rad(RA_0),DEC[src_i],ra[src_i],sigma_fwhm,ant1_theta_error,ant1_phi_error,Ant1_beam_size_error)*bs.gaussian_beam(np.deg2rad(DEC_0),np.deg2rad(RA_0),DEC[src_i],ra[src_i],sigma_fwhm,ant1_theta_error,ant1_phi_error,Ant1_beam_size_error)*np.exp(-1j*2.0*np.pi*np.dot(bl_xyz/lambda_,src_pos_lmn[src_i])) for src_i in range(len(src_pos_lmn))])
                
            vis[lst_i,:] = vis_freq
                
                
        return vis
    

def vis_sim_airydisk_beam_bl(ant1_pos,ant2_pos,dish_diam,d_block,RA_0,DEC_0,RA,DEC,src_flux_spec_index,freqs,lsts,a_theta_ant1,a_phi_ant1,theta_c_ant1,phi_c_ant1,a_theta_ant2,a_phi_ant2,theta_c_ant2,phi_c_ant2):
    
        vis =np.zeros((lsts.size,freqs.size),dtype='complex')
        
        for lst_i in range(lsts.size):
            
             #star hour angle
            ra=[lsts[lst_i] -RA[src_i] for src_i in range(RA.size)]
            
            src_pos_lmn = [pos_lmn(ra[src_i],DEC[src_i],np.deg2rad(RA_0),np.deg2rad(DEC_0)) for src_i in range(len(ra))]
            
            
            #transforming antenna position from ENU to equatorial system XYZ
            lat = np.deg2rad(DEC_0)
            xyz_1,xyz_2 = XYZ_trans_mat(ant1_pos,lat)[0],XYZ_trans_mat(ant2_pos,lat)[0]
            
            #computing the baseline
            bl_xyz = xyz_2 - xyz_1

           

            vis_freq = np.zeros(freqs.size,dtype ='complex')
            for freq_i in range(len(freqs)):
		    lambda_ = (3e8)/(freqs[freq_i]*10**6)
                    vis_freq[freq_i]= np.sum([flux_func(freqs[freq_i],src_flux_spec_index[src_i][0], src_flux_spec_index[src_i][1])*bs.param_airydisk(np.deg2rad(DEC_0),np.deg2rad(RA_0),DEC[src_i],ra[src_i],a_theta_ant1,a_phi_ant1,dish_diam,d_block,freqs[freq_i]*10**6,theta_c_ant1,phi_c_ant1)*bs.param_airydisk(np.deg2rad(DEC_0),np.deg2rad(RA_0),DEC[src_i],ra[src_i],a_theta_ant2,a_phi_ant2,dish_diam,d_block,freqs[freq_i]*10**6,theta_c_ant2,phi_c_ant2)*np.exp(-1j*2.0*np.pi*np.dot(bl_xyz/lambda_,src_pos_lmn[src_i])) for src_i in range(len(src_pos_lmn))])
                    
            vis[lst_i,:] = vis_freq
                    
                
        return vis
flux_func = lambda f,s_0,alpha : s_0*np.power(151.0/f,alpha)
       
def XYZ_trans_mat(ENU_vec,lat):
    " This matrix transform antenna position from ENU to XYZ given the latitude"
    
    trans_mat = np.matrix([[0.0,-np.sin(lat),np.cos(lat)],[1.0,0.0,0.0],[0.0,np.cos(lat),np.sin(lat)]])
    XYZ =trans_mat.dot(ENU_vec)
    return np.array(XYZ)


def baseline_proj(bl_XYZ,RA_0,Dec_0):
    "Transform baseline to projected baseline given the phase center location in RA, DEC in degrees"

    trans_matt = np.matrix([[np.sin(np.deg2rad(RA_0)),np.cos(np.deg2rad(RA_0)),0.0],[-np.sin(np.deg2rad(Dec_0))*np.cos(np.deg2rad(RA_0)),np.sin(np.deg2rad(Dec_0))*np.sin(np.deg2rad(RA_0)),np.cos(np.deg2rad(Dec_0))],[np.cos(np.deg2rad(Dec_0))*np.cos(np.deg2rad(RA_0)),-np.cos(np.deg2rad(Dec_0))*np.sin(np.deg2rad(RA_0)),np.sin(np.deg2rad(Dec_0))]])
    bl_proj = trans_matt.dot(bl_XYZ)
    
    return np.array(bl_proj)

def pos_lmn(ra,dec,ra_0,dec_0):
    
    l=np.cos(dec)*np.sin(ra-ra_0)
    m=np.sin(dec)*np.cos(dec_0)-np.cos(dec)*np.sin(dec_0)*np.cos(ra-ra_0)
    n = np.sqrt(1-(l**2 + m**2))
    
    return np.array([l,m,n])

def convert_radec_alt_az(ra,dec,telescope_lat):
    
    alt = np.arcsin(np.sin(telescope_lat)*np.sin(dec) + np.cos(telescope_lat)*np.cos(dec)*np.cos(ra))
    az  = np.arctan(np.sin(ra)/(np.cos(telescope_lat)*np.tan(dec) -np.sin(telescope_lat)*np.cos(ra)))
    return np.array([alt,az])
                    

class InterferometricArray():
    """Class that takes a list of positions and can calcuate baselines and redundancies."""
    
    def __init__(self, positions=[],ants_ind=[]):
        self.positions = np.array(positions)
        self.nant = len(positions)
        self.antNames = np.array(ants_ind)
    
    def CalculateUBLs(self, precisionFactor=100):
        """Finds the baselines, unique baselines, and related dictionaries for indexing."""
        self.blVectors, self.blNamePairs, self.blIndexPairs = [], [], []
        self.index2name = {i:name for i,name in enumerate(self.antNames)}
        self.name2index = {name:i for i,name in enumerate(self.antNames)}
        self.name2pos = {name: self.positions[self.name2index[name]] for name in self.antNames}
        for index1,ant1 in enumerate(self.antNames):
            for index2,ant2 in zip(range(index1+1,self.nant), self.antNames[index1+1:]):
                delta = np.array([int(np.round(precisionFactor*(self.positions[index1][i] - self.positions[index2][i]))) for i in range(3)])
                if delta[1] > 0 or (delta[1] == 0 and delta[0] > 0): 
                    self.blVectors.append(tuple(delta))
                    self.blNamePairs.append((ant1, ant2))
                    self.blIndexPairs.append((index1, index2))
                else: 
                    self.blVectors.append(tuple(-delta))
                    self.blNamePairs.append((ant2, ant1))
                    self.blIndexPairs.append((index2, index1))
        self.ublDict = {}
        for b in range(len(self.blVectors)):
            if self.ublDict.has_key(self.blVectors[b]): self.ublDict[self.blVectors[b]].append(self.blNamePairs[b])
            else: self.ublDict[self.blVectors[b]] = [self.blNamePairs[b]]
        self.blIndexDict = {antPair: i for i,antPair in enumerate(self.blNamePairs)}
        self.names2ublIndex = {antPair: i for i,antPairs in enumerate(self.ublDict.values()) for antPair in antPairs}
        self.indices2ublIndex = {(self.name2index[antPair[0]],self.name2index[antPair[1]]): 
                                 i for i,antPairs in enumerate(self.ublDict.values()) for antPair in antPairs}
        self.ublVectors = np.array([self.name2pos[antList[0][0]]-self.name2pos[antList[0][1]] for antList in self.ublDict.values()])
        self.ublGroups = [antList for antList in self.ublDict.values()]
        print "With", len(self.positions), "antennas there are", len(self.ublDict.items()), "unique baselines."
        self.nbl, self.nubl = len(self.blNamePairs), len(self.ublVectors)
        return [self.blNamePairs,self.blVectors,self.ublVectors,self.ublGroups,self.ublDict]
    









	
















