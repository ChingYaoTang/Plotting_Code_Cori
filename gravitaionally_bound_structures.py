import yt
from yt.data_objects.level_sets.api import *
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

# A validator function must only accept a Clump object and either return True or False.
#def _minimum_H2_fraction(clump, min_value):
#    total_H2_fraction = np.sum(clump["gas","H2_fraction"]*clump["gas","cell_mass"])/np.sum(clump["gas","cell_mass"])
#    return total_H2_fraction >= min_value

#def _max_density_threshold(clump, threshold_value):
#    max_density = np.max(clump["gas","density"].in_units("g/cm**3"))
#    return max_density >= threshold_value

# The add_validator() function adds the validator to a registry that can be accessed by the clump finder. 
#add_validator("minimum_H2_fraction", _minimum_H2_fraction)
#add_validator("max_density_threshold", _max_density_threshold)



def gravitationally_bound_structures(file_dir, file_num, output_dir, H2_minimum, Max_dens_thres, D_min_CodeUnit, D_max_CodeUnit, step, target_fields):
    print("Loading data")
    ds = yt.load("{}/DD{:0>4d}/DD{:0>4d}".format(file_dir, file_num, file_num))
    # clump finder requires a data object
    data_source = ds.all_data()
    # and a field over which the contouring is to be performed
    field = ("gas", "density")

    # data object is used to create the initial Clump object that acts as the base for clump finding
    # this one just covers the whole domain.
    master_clump = Clump(data_source, field)
    
    # validator functions can be added to determine if an individual contour should be considered a real clump
    master_clump.add_validator("gravitationally_bound", use_particles=False)
    master_clump.add_validator("min_cells", 500)
    #master_clump.add_validator("minimum_H2_fraction", H2_minimum)
    #master_clump.add_validator("max_density_threshold", ds.quan(DensCode*Max_dens_thres,"g/cm**3"))

    # set the initial minimum and maximum of the contouring field, and the step size
    c_min = DensCode*D_min_CodeUnit # DensCode = 1000 H/cc
    c_max = DensCode*D_max_CodeUnit #data_source[field].max()

    # the lower value of the contour finder will be continually multiplied by the step size
    find_clumps(master_clump, c_min, c_max, step)

    # Additional items can be added with the add_info_item() function.
    #master_clump.add_info_item("clumps_angular_momentum")
    master_clump.add_info_item("mass_weighted_jeans_mass")
    master_clump.add_info_item("total_cells")
    #master_clump.add_info_item("mass")
    master_clump.add_info_item("max_number_density")
    master_clump.add_info_item("min_number_density")
    #master_clump.add_info_item("distance_to_main_clump")

    # The master clump will represent the top of a hierarchy of clumps
    # The entire clump tree can traversed with a loop
    print("\n Clump id:")
    for clump in master_clump:
        print(clump.clump_id)
    
    # The children attribute within a Clump object contains a list of all sub-clumps.
    # Each sub-clump is also a Clump object with its own children attribute, and so on
    print("\n Clump children:")
    print(master_clump.children)
    

    # The leaves attribute of a Clump object will return a list of the individual clumps 
    # that have no children of their own (the leaf clumps)
    leaf_clumps = master_clump.leaves
    
    prj = yt.ProjectionPlot(ds, 2, field, center="c", width=(3, "pc"))
    prj.set_cmap(field=field, cmap="arbre")
    prj.annotate_clumps(leaf_clumps)
    prj.save("./clumps")

    # this will save all info items that have been calculated as well as any field values specified with the fields keyword
    fn = master_clump.save_as_dataset(filename=output_dir, fields=target_fields) 

    print("The clump information has been saved in {}.".format(fn))

    return master_clump
    

def read_clumps_info(fn):
    ds_clumps = yt.load(fn)
    
    # The tree attribute associated with the dataset provides access to the clump tree.
    print("\n All clumps id:")
    for clump in ds_clumps.tree:
        print(clump.clump_id)
    
    print("\n ds_clumps.tree.children: ")
    print( ds_clumps.tree.children )
    print("\n ds_clumps.leaves:")
    print( ds_clumps.leaves )
    
    print("\n Information of leaves:")
    for clump_obj in ds_clumps.leaves:
        clump = Clumps_Info(clump_obj)
        print("Clump {}".format(clump.clump_id))
        print("max                               : ", clump.max_density)
        print("min                               : ", clump.min_density)
        print("max density radius                : ", np.sqrt((float(clump.max_density_pos_x)-1.5)**2 + (float(clump.max_density_pos_y)-1.5)**2 + (float(clump.max_density_pos_z)-1.5)**2))
        print("avg. side length of clump         : ", clump.side_length)
        print("avg, radius of clump              : ", clump.avg_radius)
        print("total mass                        : ", clump.total_mass)
        #mass_weighted_avg_jeans_mass = np.sum(clump['grid', 'jeans_mass']*clump['grid', 'cell_mass']) / clump['grid', 'cell_mass'].sum()
        print("mass-weighted avg. jeans mass     : ", clump.jeans_mass)
        #mass_weighted_avg_dyn_time = np.sum(clump['grid', 'dynamical_time']*clump['grid', 'cell_mass']) / clump['grid', 'cell_mass'].sum()
        print("mass-weighted avg. dynamical time : ", clump.dyn_time)
        #mass_weighted_avg_temp       = np.sum(clump['grid', 'temperature']*clump['grid', 'cell_mass']) / clump['grid', 'cell_mass'].sum()
        print("mass-weighted avg. temperature    : ", clump.temp)
        #total_H2_fraction            = np.sum(clump['grid', 'H2_fraction']*clump['grid', 'cell_mass']) / clump['grid', 'cell_mass'].sum()
        print("total H2 fraction                 : ", clump.H2_fraction)
        #L_clump_x = clump['grid', 'angular_momentum_x'].sum()
        #L_clump_y = clump['grid', 'angular_momentum_y'].sum()
        #L_clump_z = clump['grid', 'angular_momentum_z'].sum()
        #L_clump   = np.sqrt( L_clump_x*L_clump_x + L_clump_y*L_clump_y + L_clump_z*L_clump_z )
        print("angular momentum                  : ", clump.L)
        print("omega                             : ", clump.omega_spin)
        print("v spin                            : ", clump.omega_spin*clump.avg_radius.in_units("cm"))
        #print("K/U                               : ", 5*clump.L**2/(6*G*(clump.total_mass.in_units("g"))**3 *clump.avg_radius.in_units("cm")) )
        print("total cell                        : ", len(clump_obj["grid","cell_volume"]))
        print("\n")

    return ds_clumps 

class Clumps_Info():
    def __init__(self, clump):
        self.clump_id    = clump.clump_id
        self.max_density = clump['grid', 'density'].max()
        self.min_density = clump['grid', 'density'].min()
        self.side_length = np.power(clump['grid', 'cell_volume'].sum().in_units("pc**3"), 1/3)
        self.avg_radius  = np.power( 3*(clump['grid', 'cell_volume'].sum().in_units("pc**3"))/(4*np.pi) , 1/3)
        self.total_mass  = clump['grid', 'cell_mass'].in_units("Msun").sum()
        self.jeans_mass  = ( np.sum(clump['grid', 'jeans_mass']*clump['grid', 'cell_mass']) / clump['grid', 'cell_mass'].sum() ).in_units("Msun")
        self.dyn_time    = ( np.sum(clump['grid', 'dynamical_time']*clump['grid', 'cell_mass']) / clump['grid', 'cell_mass'].sum() ).in_units("kyr")
        self.temp        = np.sum(clump['grid', 'temperature']*clump['grid', 'cell_mass']) / clump['grid', 'cell_mass'].sum()
        self.H2_fraction = np.sum(clump['grid', 'H2_fraction']*clump['grid', 'cell_mass']) / clump['grid', 'cell_mass'].sum()
        self.L_x         = clump['grid', 'angular_momentum_x'].sum()
        self.L_y         = clump['grid', 'angular_momentum_y'].sum()
        self.L_z         = clump['grid', 'angular_momentum_z'].sum()
        self.L           = np.sqrt( self.L_x*self.L_x + self.L_y*self.L_y + self.L_z*self.L_z )
        self.omega_spin  = 5*self.L/(2*self.total_mass.in_units("g")*self.avg_radius.in_units("cm")*self.avg_radius.in_units("cm"))
        
        self.max_density_arg   = (clump["grid","density"]==self.max_density)
        self.max_density_pos_x = clump["grid", "x"][self.max_density_arg][0].in_units("pc")
        self.max_density_pos_y = clump["grid", "y"][self.max_density_arg][0].in_units("pc") 
        self.max_density_pos_z = clump["grid", "z"][self.max_density_arg][0].in_units("pc") 

