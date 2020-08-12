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

def vis_sim_gaussian_beam_bl(ant_pos_1,ant_pos_2,dish_diam,RA_0,DEC_0,RA,DEC,flux,freqs,lsts,ant1_theta_error,ant1_phi_error, ant2_theta_error,ant2_phi_error,Ant1_beam_size_error,Ant2_beam_size_error):
    
        vis =np.zeros((lsts.size,freqs.size),dtype='complex')
        
        for lst_i in range(lsts.size):
            
            RA = [HA(lst_i,RA[i]) for i in range(len(RA))]
            sky_pos_r = [[np.cos(lalt_i)*np.cos(long_i), np.cos(lalt_i)*np.sin(long_i), np.sin(lalt_i)] for (lalt_i,long_i) in zip(DEC,RA)]
            #transforming antenna position from ENU to equatorial system XYZ
            lat = DEC_0
            RA_0 = HA(lst_i,RA_0)
            xyz_1,xyz_2 = XYZ_trans_mat(ant_pos_1,lat)[0],XYZ_trans_mat(ant_pos_2,lat)[0]
        
       
            bl_xyz = xyz_2 - xyz_1
            bl_xyz_proj = baseline_proj(bl_xyz,RA_0,DEC_0)
            #compute the visibility for all frequency
            vis_freq = np.zeros(freqs.size,dtype ='complex')
            for freq_i in range(len(freqs)):
                    
                    lambda_ = (3e8)/(freqs[freq_i]*10**6)
                    sigma_fwhm = bs.sigma_func(dish_diam,freqs[freq_i])
                    theta, phi = DEC, RA
                    
                    phase = np.array([np.exp(1j*2.0*np.pi*np.dot(bl_xyz_proj[0],sky_pos_r[r_i]/np.linalg.norm(sky_pos_r[r_i]))/lambda_)  for r_i in range(len(sky_pos_r))])
                    perceive_sky = np.array([flux[i][freq_i]*bs.gaussian_beam(theta[i],phi[i],sigma_fwhm,ant1_theta_error,ant1_phi_error,Ant1_beam_size_error)*bs.gaussian_beam(theta[i],phi[i],sigma_fwhm,ant2_theta_error,ant2_phi_error,Ant2_beam_size_error) for i in range(theta.size)])
                    vis_freq[freq_i] = np.sum(perceive_sky*phase)
                    
            vis[lst_i,:] = vis_freq
                    
                
        return vis
    

def vis_sim_airydisk_beam_bl(ant_pos_1,ant_pos_2,dish_diam,d_block,RA_0,DEC_0,RA,DEC,flux,freqs,lsts,a_theta_ant1,a_phi_ant1,theta_c_ant1,phi_c_ant1,a_theta_ant2,a_phi_ant2,theta_c_ant2,phi_c_ant2):
    
        vis =np.zeros((lsts.size,freqs.size),dtype='complex')
        
        for lst_i in range(lsts.size):
            
            RA = [HA(lst_i,RA[i]) for i in range(len(RA))]
            sky_pos_r = [[np.cos(lalt_i)*np.cos(long_i), np.cos(lalt_i)*np.sin(long_i), np.sin(lalt_i)] for (lalt_i,long_i) in zip(DEC,RA)]
            #transforming antenna position from ENU to equatorial system XYZ
            lat = DEC_0
            RA_0 = HA(lst_i,RA_0)
            xyz_1,xyz_2 = XYZ_trans_mat(ant_pos_1,lat)[0],XYZ_trans_mat(ant_pos_2,lat)[0]
        
       
            bl_xyz = xyz_2 - xyz_1
            bl_xyz_proj = baseline_proj(bl_xyz,RA_0,DEC_0)
            #compute the visibility for all frequency

            vis_freq = np.zeros(freqs.size,dtype ='complex')
            for freq_i in range(len(freqs)):
                    
                    lambda_ = (3e8)/(freqs[freq_i]*10**6)
                    sigma_fwhm = bs.sigma_func(dish_diam,freqs[freq_i])
                    theta, phi = DEC, RA
                    
                    phase = np.array([np.exp(1j*2.0*np.pi*np.dot(bl_xyz_proj[0],sky_pos_r[r_i]/np.linalg.norm(sky_pos_r[r_i]))/lambda_)  for r_i in range(len(sky_pos_r))])
                    perceive_sky = np.array([flux[i][freq_i]*bs.param_airydisk(theta[i],phi[i],a_theta_ant1,a_phi_ant1,theta_c_ant1,phi_c_ant1,dish_diam,d_block,freqs[freq_i])*bs.param_airydisk(theta[i],phi[i],a_theta_ant2,a_phi_ant2,theta_c_ant2,phi_c_ant2,dish_diam,d_block,freqs[freq_i]) for i in range(theta.size)])
                    vis_freq[freq_i] = np.sum(perceive_sky*phase)
                    
            vis[lst_i,:] = vis_freq
                    
                
        return vis
    
       
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



