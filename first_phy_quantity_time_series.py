import yt
import numpy as np
import h5py
from global_variables import *

def first_phy_quantity_time_series(ns1, ne1, nskip1, file_dir, output_dir):
    KE            = [] 
    ThE           = []
    TE            = []
    Dmax          = []
    current_t     = []
    #DivFlow      = []
    TT            = []
    H2            = []
    MachRMS       = []
    Cs            = []

    for i in range(ns1,ne1+1,nskip1):
        print("Calculating dump {0}".format(i))
        ds_ = yt.load("{}/DD{:0>4d}/DD{:0>4d}".format(file_dir,i,i))
        alld = ds_.all_data()

        KE.append(np.sum(alld["kinetic_energy"]*alld["cell_volume"]))
        ThE.append(np.sum(alld["thermal_energy"]*alld["cell_mass"]))
        TE.append(np.sum(alld["total_energy"]*alld["cell_mass"]))

        Dmax.append(np.mean(np.sort(alld["density"].flat)[-100:])) # top 100 cells

        #alldsp = alld.cut_region(["obj['spherical_radius'].in_units('pc') < 1"])
        #dens_flow_div = -np.sum((alldsp["density"]*alldsp["velocity_divergence"]
        #                        +alldsp["velocity_x"]*alldsp["density_gradient_x"]
        #                        +alldsp["velocity_y"]*alldsp["density_gradient_y"]
        #                        +alldsp["velocity_z"]*alldsp["density_gradient_z"])*alldsp["cell_volume"])
        #DivFlow.append(dens_flow_div/SunMass*yrs)

        #alld1e4 = alld.cut_region(["(obj['number_density'] < 1.25e4) & (obj['number_density'] > 0.75e4)"])
        TT.append(np.sum(alld["temperature"]*alld["density"])/np.sum(alld["density"]))
        H2.append(np.sum(alld["H2_fraction"]*alld["cell_mass"])/np.sum(alld["cell_mass"]))
        MachRMS.append(np.sqrt(np.mean(np.power(alld["velocity_magnitude"]/alld['sound_speed'],2))))
        Cs.append(np.sqrt(np.mean(np.power(alld['sound_speed'],2))))

        current_t.append(ds_.current_time.in_units('Myr'))

    f = h5py.File(output_dir, mode="w")
    f.create_dataset("/KE", data=np.array(KE))
    f.create_dataset("/ThE", data=np.array(ThE))
    f.create_dataset("/TE", data=np.array(TE))
    f.create_dataset("/Dmax_arr", data=np.array(Dmax))
    #f.create_dataset("/DivFlow", data=np.array(DivFlow))
    f.create_dataset("/TT", data=np.array(TT))
    f.create_dataset("/H2", data=np.array(H2))
    f.create_dataset("/MachRMS", data=np.array(MachRMS))
    f.create_dataset("/Cs", data=np.array(Cs))

    f.create_dataset("/current_t_arr", data=np.array(current_t))
    f.close()


def max_density_time_series(ns1, ne1, nskip1, file_dir, output_dir):
    Dmax          = []
    current_t     = []

    for i in range(ns1,ne1+1,nskip1):
        print("Calculating dump {0}".format(i))
        ds_ = yt.load("{}/DD{:0>4d}/DD{:0>4d}".format(file_dir,i,i))
        alld = ds_.all_data()

        Dmax.append(np.mean(np.sort(alld["density"].flat)[-10:])) # top 100 cells
        current_t.append(ds_.current_time.in_units('Myr'))

    f = h5py.File(output_dir, mode="w")
    f.create_dataset("/Dmax_arr", data=np.array(Dmax))
    f.create_dataset("/current_t_arr", data=np.array(current_t))
    f.close()
