import yt
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib as mplb
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy import optimize as op
from scipy import stats
import h5py
import os

from gravitaionally_bound_structures import *

def find_the_field_extrema(file_dir, field):
    ds = yt.load( file_dir )
    min_field_z, max_field_z = ds.all_data().quantities.extrema(field)
    return [min_field_z, max_field_z]
    
##########################################################################################################################

def phase_diagram(Files, output_dir, ncol, nrow, figsize_x, figsize_y, label_size, tick_size, model_name,
                  field_x, field_x_label, 
                  field_y, field_y_label, 
                  field_z, field_z_label, 
                  #min_field_z, max_field_z,
                  weight_field=None, colormap='jet'):

    if (len(Files) != ncol*nrow):
        print("Length of Files does not equal to the total number of plots.\n")
        return 0
    
    # Load all the data
    ds_list = []
    for file in Files:
        ds_list.append(yt.load( file ))

    # Find the upper and lower limits of each field
    min_field_x, max_field_x = ds_list[0].all_data().quantities.extrema(field_x)
    min_field_y, max_field_y = ds_list[0].all_data().quantities.extrema(field_y)
    min_field_z, max_field_z = ds_list[0].all_data().quantities.extrema(field_z)
    for i in range(1, len(Files)):
        min_field_x_tmp, max_field_x_tmp = ds_list[i].all_data().quantities.extrema(field_x)
        min_field_x          = np.min([min_field_x_tmp, min_field_x])
        max_field_x          = np.max([max_field_x_tmp, max_field_x])

        min_field_y_tmp, max_field_y_tmp = ds_list[i].all_data().quantities.extrema(field_y)
        min_field_y          = np.min([min_field_y_tmp, min_field_y])
        max_field_y          = np.max([max_field_y_tmp, max_field_y])

        min_field_z_tmp, max_field_z_tmp = ds_list[i].all_data().quantities.extrema(field_z)
        min_field_z          = np.min([min_field_z_tmp, min_field_z])
        max_field_z          = np.max([max_field_z_tmp, max_field_z])
        
    # Create a blank figure
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
    grid = AxesGrid(fig, rect=(0.1,0.1,0.8,0.8),
                nrows_ncols = (ncol, nrow),
                axes_pad = 0.35,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="2%",
                aspect=False
                )

    for i, ds in enumerate(ds_list):
        # Load the data and create a single plot
        ad = ds.all_data()
        p  = yt.PhasePlot( ad, field_x, field_y, [field_z], weight_field=weight_field)
    
        # Ensure the axes and colorbar limits match for all plots
        p.set_xlim( min_field_x, max_field_x )
        p.set_ylim( min_field_y, max_field_y )
        p.set_zlim( field_z, min_field_z, max_field_z )
        #p.annotate_title('({}) Model {}'.format(chr(97+i),case_ds[i]))
        p.set_cmap(field=field_z, cmap=colormap)
        #p.set_figure_size(4)

        # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = p.plots[(field_z)]
        plot.figure = fig
        plot.axes = grid[i].axes
        if i == 0:
            plot.cax = grid.cbar_axes[i]
        
        #grid[i].axes.yaxis.set_visible(False)
        #grid[i].axes.xaxis.set_visible(False)
        
        # Actually redraws the plot.
        p._setup_plots()

        # Modify the axes properties **after** p._setup_plots() so that they are not overwritten.
        #plot.axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]) )
   
    for i in range(ncol*nrow):
        grid.cbar_axes[i].set_ylabel(field_z_label, fontsize=label_size)
        grid.cbar_axes[i].tick_params(labelsize = label_size, color='k', direction='out', which='both', length=5, width=1)
        #cb.ax.minorticks_on()
        
        ### model name and time of the snapshot
        # grid[i].axes.text(-0.45, 0.4, alphabat[i] + ' ' + model_name[i], fontsize=text_size, color='w')

        grid[i].axes.set_ylabel(field_y_label, fontsize=label_size)
        grid[i].tick_params(axis="y", labelsize=tick_size, direction='in', which='both', length=5, width=1)
        grid[i].tick_params(axis="y", which='minor', length=0, width=0)

        grid[i].axes.set_xlabel(field_x_label, fontsize=label_size)
        grid[i].tick_params(axis="x", labelsize=tick_size, direction='in', which='both', length=5, width=1 )
        '''
        if i%ncol == 0:
            #grid[i].axes.yaxis.set_visible(True)
            #grid[i].axes.yaxis.tick_left()
            grid[i].tick_params(axis="y", labelsize=ylabel_size, direction='out', which='both')
            grid[i].tick_params(axis="y", which='minor', length=7, width=2)
        if i >= (nrow-1)*ncol:
            print("xx")
            #grid[i].axes.xaxis.set_visible(True)
            #grid[i].axes.xaxis.tick_bottom()
            grid[i].tick_params(axis="x", labelsize=xlabel_size, direction='out', which='both')
            grid[i].tick_params(axis="x", which='minor', length=7, width=2)
        '''
        grid[i].axes.set_title( "({}) {}".format( chr(97+i), model_name[i]), fontsize=label_size, color='k' )
    
    plt.savefig( output_dir )

########################################################################################################################

def clump_phase_diagram(Clumps, Files, output_dir, ncol, nrow, figsize_x, figsize_y, label_size, tick_size, model_name,
                        field_x, field_x_label,
                        field_y, field_y_label,
                        field_z, field_z_label,
                        title_size, xticks, yticks,
                        weight_field=None, colormap='jet'):

    if (len(Clumps) != ncol*nrow):
        print("Length of Files does not equal to the total number of plots.\n")
        return 0


    # Find the upper and lower limits of each field
    min_field_x = Clumps[0]["grid", field_x].min()
    max_field_x = Clumps[0]["grid", field_x].max()
    min_field_y = Clumps[0]["grid", field_y].min()
    max_field_y = Clumps[0]["grid", field_y].max()
    min_field_z = Clumps[0]["grid", field_z].min()
    max_field_z = Clumps[0]["grid", field_z].max()
    for i in range(1, len(Clumps)):
        min_field_x_tmp      = Clumps[i]["grid", field_x].min()
        max_field_x_tmp      = Clumps[i]["grid", field_x].max()
        min_field_x          = np.min([min_field_x_tmp, min_field_x])
        max_field_x          = np.max([max_field_x_tmp, max_field_x])

        min_field_y_tmp      = Clumps[i]["grid", field_y].min()
        max_field_y_tmp      = Clumps[i]["grid", field_y].max()
        min_field_y          = np.min([min_field_y_tmp, min_field_y])
        max_field_y          = np.max([max_field_y_tmp, max_field_y])

        min_field_z_tmp      = Clumps[i]["grid", field_z].min()
        max_field_z_tmp      = Clumps[i]["grid", field_z].max()
        min_field_z          = np.min([min_field_z_tmp, min_field_z])
        max_field_z          = np.max([max_field_z_tmp, max_field_z])

    # Create a blank figure
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
    grid = AxesGrid(fig, rect=(0.08,0.08,0.8,0.88),
                    nrows_ncols = (nrow, ncol),
                    axes_pad = 0.4,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="3%",
                    cbar_pad="2%",
                    aspect = False
                   )

    ds_list = []
    for file in Files:
        ds_list.append(yt.load( file ))

    for i, ds in enumerate(ds_list):
        clump = Clumps_Info(Clumps[i])
        print("clump max density position = ", [clump.max_density_pos_x.in_units("cm"),
                                                clump.max_density_pos_y.in_units("cm"),
                                                clump.max_density_pos_z.in_units("cm")])
        # Load the data and create a single plot
        sp = ds.sphere((float(clump.max_density_pos_x.in_units("code_length")),
                        float(clump.max_density_pos_y.in_units("code_length")),
                        float(clump.max_density_pos_z.in_units("code_length"))),
                       (float(clump.avg_radius.in_units("pc"))*0.9, "pc"))
        print("sphere center position     = ",sp)
        p  = yt.PhasePlot( sp, field_x, field_y, [field_z], weight_field=weight_field)

        # Ensure the axes and colorbar limits match for all plots
        p.set_xlim( min_field_x, max_field_x )
        p.set_ylim( min_field_y, max_field_y )
        p.set_zlim( field_z, min_field_z, max_field_z )
        p.set_cmap(field=field_z, cmap=colormap)
        p.set_log(field_y, False)
        #p.set_log(field_z, False)
        #p.set_figure_size(4)

        # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = p.plots[(field_z)]
        plot.figure = fig
        #grid[i].cax.minorticks_on()
        plot.axes = grid[i].axes
        if i == 0:
            plot.cax = grid.cbar_axes[i]

        #grid[i].axes.yaxis.set_visible(False)
        #grid[i].axes.xaxis.set_visible(False)

        # Actually redraws the plot.
        p._setup_plots()

        # Modify the axes properties **after** p._setup_plots() so that they are not overwritten.
        #plot.axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]) )


        grid.cbar_axes[i].set_ylabel(field_z_label, fontsize=label_size)
        grid.cbar_axes[i].tick_params(labelsize = tick_size, color='k', direction='out', which='both', length=5, width=1)
        #cb.ax.minorticks_on()

        ### model name and time of the snapshot
        # grid[i].axes.text(-0.45, 0.4, alphabat[i] + ' ' + model_name[i], fontsize=text_size, color='w')
        #grid[i].minorticks_on()
        grid[i].axes.yaxis.tick_left()
        grid[i].axes.set_ylabel(" ", fontsize=label_size)
        grid[i].set_yticks(yticks)
        grid[i].tick_params(axis="y", labelsize=tick_size, direction='out', which='both', length=5, width=1.2)
        grid[i].tick_params(axis="y", which='minor', length=3, width=0.8)
        #grid[i].yaxis.set_minor_locator( mplb.ticker.AutoMinorLocator(10) )

        grid[i].axes.xaxis.tick_bottom()
        grid[i].axes.set_xlabel(" ", fontsize=label_size)
        grid[i].set_xticks(xticks)
        grid[i].tick_params(axis="x", labelsize=tick_size, direction='out', which='both', length=5, width=1.2 )
        grid[i].tick_params(axis="x", which='minor', length=3, width=0.8)
        #grid[i].xaxis.set_minor_locator( mplb.ticker.AutoMinorLocator(10) )

        #grid[i].axes.set_title( "{}".format(model_name[i]), fontsize=title_size, color='k' )
        grid[i].axes.set_title( "({}) {}".format( chr(97+i), model_name[i]), fontsize=title_size, color='k' )

    fig.text(0.5, 0.0075, field_x_label, ha='center', fontsize=label_size)
    fig.text(0.005, 0.5, field_y_label, va='center', rotation='vertical', fontsize=label_size)
    plt.savefig( output_dir, dpi=300 )

