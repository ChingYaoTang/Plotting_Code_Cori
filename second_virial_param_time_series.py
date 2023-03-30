import yt
import numpy as np
import h5py
import time
from enum import Enum
from global_variables import *


# Calculating enclosed mass of NFW halo
def M_NFW_enc(r, M_vir, R_vir, c):
    return M_vir * (np.log(1 + r*c/R_vir) - (r*c/R_vir)/(1 + r*c/R_vir)) / (np.log(1 + c) - (c)/(1 + c))

def second_virial_param_time_series(ns2, ne2, nskip2, ThisHalo, file_dir, output_dir):
    KE              = []
    ThE             = []
    V_NFW           = []
    V_Self_Shell    = []
    V_Self_Green    = []
    Es              = []
    Beta            = []
    Menc            = []
    #Vel_k           = []
    #Vel_k2          = []
    VelRMS          = []
    MachRMS       = []
    current_t       = []
    #k0_             = 4

    #green_r        = np.array([[[1/np.sqrt(i*i+j*j+k*k) for i in range(rsl)] for j in range(rsl)] for k in range(rsl)])
    """
    g_extend        = np.array([[[1/np.sqrt( np.power((rsl-np.abs(i-rsl)),2) + 
                                             np.power((rsl-np.abs(j-rsl)),2) + 
                                             np.power((rsl-np.abs(k-rsl)),2) )
                                   for i in range(2*rsl)]
                                 for j in range(2*rsl)]
                               for k in range(2*rsl)] )
    """
    start_time      = time.time()
    g_extend_i      = np.array([[[ ( i )
                                   for i in range(2*rsl)]
                                 for j in range(2*rsl)]
                               for k in range(2*rsl)] )
    mid_time        = time.time()
    print("g_extend_i   costs {:>10.6f} sec to be done.".format(mid_time-start_time))
    g_extend        = 1/np.sqrt( np.power( (rsl-np.abs(g_extend_i                  - rsl)),2) + 
                                 np.power( (rsl-np.abs(g_extend_i.transpose(0,2,1) - rsl)),2) + 
                                 np.power( (rsl-np.abs(g_extend_i.transpose(2,1,0) - rsl)),2) )
    end_time        = time.time()
    print("g_extend     costs {:>10.6f} sec to be done.".format(end_time-start_time))


    g_extend[0,0,0] = 3
    #g_extend       /= 3.61602890625e+16
    start_time      = time.time()
    g_FT            = np.fft.rfftn(g_extend)
    end_time        = time.time()
    print("g_FT         costs {:>10.6f} sec to be done.\n".format(end_time-start_time))

    # First datadump has no AMR structure yet, and the index of ['spherical_radius'] is correct
    datadump_0_dir               = file_dir.split('/data')[0]
    ds_                          = yt.load("{}/data/DD{:0>4d}/DD{:0>4d}".format(datadump_0_dir,0,0))
    alld                         = ds_.all_data()
    # The data array of all_data() is 1D, so we need to reshape it 
    spherical_system_idx         = np.array(alld['index','spherical_radius'].in_units('cm') < Rc).reshape(256,256,256)
    spherical_radius_arr         = np.array(alld['index','spherical_radius'].in_units('cm')).reshape(256,256,256)
    cell_width                   = float(alld['index','dx'][0].in_units('cm'))

    for i in range(ns2,ne2+1,nskip2):
        print("Calculating dump {0}".format(i))
        # Timing for 1 loop
        start_time               = time.time()
        ds_                      = yt.load("{}/DD{:0>4d}/DD{:0>4d}".format(file_dir,i,i))
        #alld                      = ds_.all_data()
        alld                     = ds_.covering_grid(level=0, left_edge=[0.0,0.0,0.0], dims=ds_.domain_dimensions)
        # Spherical system with radius Rc
        #spherical_system_idx     = (alld['index','spherical_radius'].in_units('cm') < Rc)
        #sphd = alld.cut_region(["obj['spherical_radius'].in_units('cm') < {0}".format(Rc)])
    
        # Kinetic Energy within Rc
        KE.append(float(np.sum( (alld["kinetic_energy"]*alld["cell_volume"])[spherical_system_idx] )))
    
        # Thermal Energy within Rc
        ThE.append(float(np.sum( (alld["thermal_energy"]*alld["cell_mass"])[spherical_system_idx] )))
    
        # Gravitational Potential Energy within Rc
        # Spherical symmetric calculation
        cell_mass_arr            = np.array(alld['gas','cell_mass'].in_units('g'))
    
        r                        = cell_width
        sphere_idx               = spherical_radius_arr < r
        Menc_self                = (np.sum(cell_mass_arr[sphere_idx]))
        # r' < cell_width
        potential_self           = -3*G                                                                *Menc_self**2/(5*r)
        potential_NFW            = -G*M_NFW_enc(r, ThisHalo.value[0], ThisHalo.value[1], concentration)*Menc_self/r

        while r < Rc:
            r                   += cell_width
            #print("r =  {0}, cell width = {1}".format(r,cell_width))
            sphere_idx           = (
                                    (spherical_radius_arr >= r-cell_width) & 
                                    (spherical_radius_arr < r))
            # Shell mass within [ r-cell_width, r )
            Mshell               = (np.sum(cell_mass_arr[sphere_idx]))
            potential_self      += -G*Menc_self                                                                   *Mshell/(r-0.5*cell_width)
            potential_NFW       += -G*M_NFW_enc(r-cell_width, ThisHalo.value[0], ThisHalo.value[1], concentration)*Mshell/(r-0.5*cell_width)
            Menc_self           += Mshell
        V_Self_Shell.append(potential_self)
        V_NFW.append(potential_NFW)
        
        # Calculating potential by Green function
        M_extend                 = np.zeros((512,512,512))
        M_extend[:256,:256,:256] = cell_mass_arr
        # Timimg for FT
        mid_time                 = time.time()
        M_FT                     = np.fft.rfftn(M_extend)
        # Only Phi[:256,:256,:256] (quadrant 2) is the ture potential values
        Phi                      = -G/cell_width* (np.fft.irfftn(g_FT*M_FT))[:256,:256,:256]
        # Timing fot FT
        mid1_time                = time.time()
        print("Phi and M_FT cost  {:>10.6f} sec to be done.".format(mid1_time-mid_time))
        # Calculate the potential within the shperical system
        V_Self_Green.append( 0.5 * np.sum(Phi[spherical_system_idx]
                                      * cell_mass_arr[spherical_system_idx]) )
    
    
        # Surface Pressure
        sphere_idx               = (
                                    (spherical_radius_arr >= Rc-cell_width) & 
                                    (spherical_radius_arr < Rc) )
        Pavg                     = float(np.average(alld['gas','pressure'][sphere_idx]))
        Es.append(Pavg*4*np.pi*np.power(Rc,3))

        # Beta
        Beta.append((3*(gamma-1)*ThE[-1]+2*KE[-1])/(Es[-1]-V_NFW[-1]-V_Self_Green[-1])-1)
        
        # Menc
        Menc.append(np.sum(cell_mass_arr[spherical_system_idx]))
        
        # Vel
        Vel_arr                  = np.array(alld["velocity_magnitude"].in_units('cm/s'))
        VelRMS.append( np.sqrt(np.mean(np.power(Vel_arr,2))) )
        #Vel_FT                   = np.fft.rfftn(Vel_arr)
        #Vel_k.append( (Vel_FT[k0_,0,0]+Vel_FT[0,k0_,0]+Vel_FT[0,0,k0_])/3 )
        #Vel_k2.append( (Vel_FT[2,0,0]+Vel_FT[0,2,0]+Vel_FT[0,0,2])/3 )
        MachRMS.append(np.sqrt(np.mean(np.power(alld["velocity_magnitude"]/alld['sound_speed'],2))))
        
        # Time
        current_t.append(ds_.current_time.in_units('Myr'))
    
        # Timing for 1 loop
        end_time   = time.time()
        print("Dump {:>3d}     costs {:>10.6f} sec to be done. \n".format(i,end_time-start_time))
        
    f = h5py.File(output_dir, mode="w")
    f.create_dataset("/KE", data=np.array(KE))
    f.create_dataset("/ThE", data=np.array(ThE))
    f.create_dataset("/V_NFW", data=np.array(V_NFW))
    f.create_dataset("/V_Self_Shell", data=np.array(V_Self_Shell))
    f.create_dataset("/V_Self_Green", data=np.array(V_Self_Green))
    f.create_dataset("/Es", data=np.array(Es))
    f.create_dataset("/Beta", data=np.array(Beta))
    f.create_dataset("/Menc", data=np.array(Menc))
    #f.create_dataset("/Vel_k", data=np.array(Vel_k))
    #f.create_dataset("/Vel_k2", data=np.array(Vel_k2))
    f.create_dataset("/VelRMS", data=np.array(VelRMS))
    f.create_dataset("/MachRMS", data=np.array(MachRMS))
    f.create_dataset("/current_t_arr", data=np.array(current_t))
    f.close()
