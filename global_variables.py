import yt
from enum import Enum

## Basical units and physical parameters
rsl                  = 256
ds_units             = yt.load("/global/cscratch1/sd/cytang/testM2C2S/data/DD{:0>4d}/DD{:0>4d}".format(0,0))
LenCode              = float(ds_units.length_unit.in_units("cm"))          # 3 pc
DensCode             = float(ds_units.mass_unit.in_units("g"))/LenCode**3  # 1000 H/cc
TimeCode              = float(ds_units.time_unit.in_units("Myr"))
pc_cm                = float(ds_units.quan(1, "pc").in_units("cm"))
Msun_g               = float(ds_units.quan(1, "Msun").in_units("g"))
yr_s                 = float(ds_units.quan(1, "yr").in_units("s"))
G                    = 6.67408e-8                                          # G constant[cm3 g-1 s-2]
gamma                = 1.6667                                              # Adiabatic index
kb                   = 1.3806504e-16                                       # Boltzmann's constant [cm2 g s-2 K-1] or [erg K-1]
mu                   = 1.22                                                # Average
mh                   = 1.67262171e-24                                      # Mass of hydrogen [g]
del ds_units

## Radius of cloud system
Rc_pc                = 1.5
Rc                   = Rc_pc  * pc_cm     # spherical system radius
#Rc                   = 2.587927 * pc_cm   # maximum distance from the center

## Stochastic forcing parameters
alpha                = 2

## NFW parameters
concentration        = 20
class MiniHalo(Enum):
    Open    = [3e5   * Msun_g , 110 * pc_cm]
    Compact = [3.1e5 * Msun_g , 97  * pc_cm]

## yt filters
pop3_type = 5
def stars(pfilter, data):
    filter = data[(pfilter.filtered_type, "particle_type")] == pop3_type
    return filter


