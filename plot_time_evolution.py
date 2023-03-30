import yt
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib as mplb
from scipy import optimize as op
from scipy import stats
import h5py
import os

from global_variables import *

def plot_Dmax_time(file_dir, output_dir):
    f = h5py.File(file_dir, mode="r")

    plt.figure(figsize=(10,10))
    plt.plot(f['current_t_arr'][:],f['Dmax_arr'][:], 'r-')
    plt.title("")
    plt.grid()
    #plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='both',useMathText=True,useOffset=True)
    plt.minorticks_on()
    plt.xlabel("Time (Myr)",size=16)
    plt.ylabel(r"$\rho_{max} \, (\frac{g}{cm^3})$",size=16)

    plt.savefig(output_dir)
    plt.show()

    f.close()

def plot_Tavg_time(file_dir, output_dir, nskip1):
    f = h5py.File(file_dir, mode="r")

    plt.figure(figsize=(10,10))
    plt.plot(f['current_t_arr'][:],f['TT'][:], 'r-')
    plt.title("")
    plt.grid()
    #plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='both',useMathText=True,useOffset=True)
    plt.minorticks_on()
    plt.xlabel("Time (Myr)",size=16)
    plt.ylabel("T (K)",size=16)

    plt.savefig(output_dir)
    plt.show()

    for i in range(10):
        print("Average T = {:.7f}, Physical time = {:.7f}, Dump number = {:2d}".format(
              (np.sort(f['TT'][1:]))[i], f['current_t_arr'][1:][np.argsort(f['TT'][1:])][i], (np.argsort(f['TT'][1:])+1)[i]*nskip1))

    f.close()

def plot_H2_time(file_dir, output_dir):
    f = h5py.File(file_dir, mode="r")

    plt.figure(figsize=(10,10))
    plt.plot(f['current_t_arr'][:],f['H2'][:], 'b-')
    plt.title("")
    plt.grid()
    #plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='both',useMathText=True,useOffset=True)
    plt.minorticks_on()
    plt.xlabel("Time (Myr)",size=16)
    plt.ylabel(r"$f_{H_2}$",size=16)
    plt.tick_params(axis='both',labelsize=16)

    plt.savefig(output_dir)
    plt.show()

    f.close()

def plot_MachRMS_time(file_dir, output_dir):
    f = h5py.File(file_dir, mode="r")

    plt.figure(figsize=(10,10))
    plt.plot(f['current_t_arr'][:],f['MachRMS'][:], 'k-')
    plt.title("")
    plt.grid()
    #plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='both',useMathText=True,useOffset=True)
    plt.minorticks_on()
    plt.xlabel("Time (Myr)",size=16)
    plt.ylabel(r"$\mathcal{M}_{rms}$",size=16)
    plt.tick_params(axis='both',labelsize=16)

    plt.savefig(output_dir)
    plt.show()

    f.close()

def plot_Cs_time(file_dir, output_dir):
    f = h5py.File(file_dir, mode="r")

    plt.figure(figsize=(10,10))
    plt.plot(f['current_t_arr'][:],f['Cs'][:]/1e5, 'g-')
    plt.title("")
    plt.grid()
    #plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='both',useMathText=True,useOffset=True)
    plt.minorticks_on()
    plt.xlabel("Time (Myr)",size=16)
    plt.ylabel(r"$\mathcal{C}_{S_{rms}} (\frac{km}{s})$",size=16)
    plt.tick_params(axis='both',labelsize=16)
    
    plt.savefig(output_dir)
    plt.show()

    f.close()

def show_target_Beta(file_dir, critical, nskip2, ns2):
    f = h5py.File(file_dir, mode="r")
    
    idx_0 = 0
    while f['Beta'][idx_0] > critical and idx_0 < len(f['Beta'][:])-1:
        idx_0 += 1
        #print(idx_0)
    for i in range(  np.min( [9, len(f['Beta'][:])-idx_0] ) +8  ):
        print("Beta = {:>10.7f}, MachRMS = {:>6.3f}, VelRMS = {:>6.3e}, Physical time = {:>6.3f}, Dump number = {:>2d}".format(
            f['Beta'][idx_0-8+i],
            f['MachRMS'][idx_0-8+i], f['VelRMS'][idx_0-8+i],
            f['current_t_arr'][idx_0-8+i], (idx_0-8+i)*nskip2+ns2) 
        )
    print("\n")
    
    f.close()

def plot_Beta_Menc_time(file_dir, output_dir, nskip2, ns2):
    f = h5py.File(file_dir, mode="r")
    
    fig, ax1 = plt.subplots(figsize=(10,10),sharex=True,sharey=True)
    ax2 = ax1.twinx()

    plot1  = ax1.plot(f['current_t_arr'][:], f['Beta'][:], 'r-',label=r"$\beta$")
    ax1.set_ylabel(r"$\beta$", color='red',size=16)
    ax1.set_xlabel("Time (Myr)",size=16) 
    ax1.tick_params(axis='y', labelcolor='red')
    ax1.grid()
    ax1.minorticks_on()

    plot2 = ax2.plot(f['current_t_arr'][:], f['Menc'][:]/Msun_g, 'b-', alpha=1, label=r"$M_{enc}$")
    ax2.set_ylabel(r"$M_{enc} (M_\odot)$", color='blue',size=16)
    ax2.set_xlabel("Time (Myr)",size=16)
    ax2.tick_params(axis='y', labelcolor='blue')
    ax2.grid()
    ax2.minorticks_on()

    lines = plot1+plot2
    ax1.legend(lines, [l.get_label() for l in lines], loc="center right",fontsize = 12)

    plt.savefig(output_dir)
    plt.show()

    """
    for i in range( np.min( [5, len(f['Beta'][:])] ) ):
        print("Beta = {:>10.7f}, Physical time = {:>10.7f}, Dump number = {:>2d}".format(
            (np.sort(f['Beta'][:]))[i], f['current_t_arr'][:][np.argsort(f['Beta'][:])][i], (np.argsort(f['Beta'][:]))[i]*nskip2+ns2))
    """
    idx_0 = 0
    while f['Beta'][idx_0] > 0 and idx_0 < len(f['Beta'][:])-1:
        idx_0 += 1
        #print(idx_0)
    for i in range(  np.min( [4, len(f['Beta'][:])-idx_0] )+1  ):
        print("Beta = {:>10.7f}, Physical time = {:>10.7f}, Dump number = {:>2d}".format(
            f['Beta'][idx_0-1+i], f['current_t_arr'][idx_0-1+i], (idx_0-1+i)*nskip2+ns2) )
    print("\n")
    f.close()


def plot_Tturnover_VelRMS_time(file_dir, output_dir):
    f = h5py.File(file_dir, mode="r")
  
    T_turnover = (LenCode/alpha) / f['VelRMS'][:]/yr_s/1e6
    fig, ax1 = plt.subplots(figsize=(10,10),sharex=True,sharey=True)
    ax2 = ax1.twinx()
  
    plot1  = ax1.plot(f['current_t_arr'][:], T_turnover, 'b-',label=r"$T_{turnover}$")
    ax1.set_ylabel(r"$T_{turnover} (Myrs)$", color='blue',size=16)
    ax1.set_xlabel("Time (Myr)",size=16) 
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid()
    ax1.minorticks_on()

    plot2 = ax2.plot(f['current_t_arr'][:], f['VelRMS'][:]/1e5, 'k-', alpha=1, label=r"$v_{RMS}$")
    ax2.set_ylabel(r"$v_{RMS} (km/s)$", color='black',size=16)
    ax2.set_xlabel("Time (Myr)",size=16)
    ax2.tick_params(axis='y', labelcolor='black')
    ax2.grid()
    ax2.minorticks_on()

    lines = plot1+plot2
    ax1.legend(lines, [l.get_label() for l in lines], loc="center right",fontsize = 12)

    plt.savefig(output_dir)
    plt.show()

    i        = 1
    dumpnum  = int(round( f['current_t_arr'][i]/(f['current_t_arr'][1] - f['current_t_arr'][0]) ))
    #print(f['current_t_arr'][i]/(f['current_t_arr'][1] - f['current_t_arr'][0]))
    duration = f['current_t_arr'][:] - f['current_t_arr'][0]
    OneTurnover = 0
    try:
        while duration[i] < 2*T_turnover[i]:
            print("At dump {:04d}, duration = {:<5.3f}, turnover time = {:<5.3f}"
                .format(dumpnum, duration[i], T_turnover[i]))
            if duration[i] > T_turnover[i] and OneTurnover == 0:
                OneTurnover = 1
                print("\nAt {:<5.3f} Myrs (dump {:04d}), the decaying turbulence has experienced {:<5.3f} Myrs (delta dump ={: 3d}), which is greater than a turnover time {:<5.3f} Myrs.\n"
                    .format(f['current_t_arr'][i],
                        int(np.ceil(dumpnum)),
                        duration[i],
                        i,
                        T_turnover[i]) )
            dumpnum += 1
            i       += 1

        print("\nAt {:<5.3f} Myrs (dump {:04d}), the decaying turbulence has experienced {:<5.3f} Myrs (delta dump ={: 3d}), which is greater than twice the turnover time {:<5.3f} Myrs.\n"
            .format(f['current_t_arr'][i],
                    int(np.ceil(dumpnum)), 
                    duration[i], 
                    i,
                    T_turnover[i]) )
    except:
        print("\nThe duration of total dumps is not enough to exceed twice the turnover time.\n")
    f.close()

def plot_potentials_time(file_dir, output_dir):
    f = h5py.File(file_dir, mode="r")

    fig, ax1 = plt.subplots(figsize=(10,10),sharex=True,sharey=True)
    ax2 = ax1.twinx()

    plot  = ax1.plot(f['current_t_arr'][:],f['V_NFW'][:]/Msun_g, 'b-',label=r"$U_{NFW}$")
    plot1 = ax1.plot(f['current_t_arr'][:],f['V_Self_Shell'][:]/Msun_g, 'b--',label=r"$U_{self,sphere}$")
    plot2 = ax1.plot(f['current_t_arr'][:],f['V_Self_Green'][:]/Msun_g, 'b-.',label=r"$U_{self,FT}$")
    ax1.set_ylabel(r"$U (erg)$", color='blue',size=16)
    ax1.set_xlabel("Time (Myr)",size=16) 
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid()
    ax1.minorticks_on()

    plot3 = ax2.plot(f['current_t_arr'][:],(f['V_Self_Shell'][:]-f['V_Self_Green'][:])/f['V_Self_Green'][:]*100, 
                        'k-', alpha=1, label=r"Error of $U_{self}$")
    ax2.set_ylabel('Error (%)', color='black',size=16)
    ax2.set_xlabel("Time (Myr)",size=16)
    ax2.tick_params(axis='y', labelcolor='black')
    ax2.grid()
    ax2.minorticks_on()

    lines = plot+plot1+plot2+plot3
    ax1.legend(lines, [l.get_label() for l in lines], loc="lower center",fontsize = 12)

    plt.savefig(output_dir)
    plt.show()

    f.close()


def slice_plot(file_dir, output_dir, ns, ne, nskip, field_name, cmap_name, annotate, enclosed_mass=0):
    for i in range(ns, ne+1, nskip):                          
        ds = yt.load("{}/DD{:0>4d}/DD{:0>4d}".format(file_dir,i,i))
        slc = yt.SlicePlot(ds, "z", field_name, center="c")
        slc.set_cmap(field=field_name, cmap=cmap_name)
        slc.annotate_timestamp(corner='upper_left', draw_inset_box=True)
        
        if annotate == "velocity":
            slc.annotate_velocity(plot_args={"headwidth": 3,"color":"magenta"})
        elif annotate == "streamlines":
            slc.annotate_streamlines(("gas", "velocity_x"), ("gas", "velocity_y"))
        
        if enclosed_mass == 1:
            sp          = ds.sphere("c", (1.0, "pc"))
            baryon_mass = sp.quantities.total_quantity( ("gas", "cell_mass") )
            slc.annotate_text([-1.1,-1.4], "Total gas mass within 1pc sphere = %7.3f $M_\odot$" % baryon_mass.in_units("Msun"), coord_system='plot')
        

        slc.annotate_title('{} slice plot'.format(field_name))
        ad = ds.all_data()

        max_dens     = ad.quantities.extrema("density")[1]
        rad_max_dens = ad["spherical_radius"][ad["density"]==max_dens][0]
        fH2_max_dens = ad["H2_p0_fraction"][ad["density"]==max_dens][0]
        T_max_dens   = ad["temperature"][ad["density"]==max_dens][0]

        slc.annotate_text([-1.1,-1.0], "radius of max density = %4.3f pc" % rad_max_dens.in_units('pc'), coord_system='plot')
        slc.annotate_text([-1.1,-1.1], "max density = %0.6e g/$cm^3$" % max_dens.in_units('g/cm**3'), coord_system='plot')
        slc.annotate_text([-1.1,-1.2], "fH2 of max density = %0.6e" % fH2_max_dens, coord_system='plot')
        slc.annotate_text([-1.1,-1.3], "T of max density = %7.3f K" % T_max_dens, coord_system='plot')
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        slc.save("{}/{}{}_{:0>4d}".format(output_dir, field_name, file_dir.split('/data')[-1], i))



