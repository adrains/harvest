"""Script for vetting TOIs using Gaia DR2 isochones 
"""
from __future__ import division, print_function
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from scipy.interpolate import interp1d

isochrones_file = "padova_isochrone.csv"
isochrones_file = "padova_isochrone_age.dat"
#sochrones_file = "padova_isochrone_metallicity.dat"
toi_file = "tois.csv"

names = ["Zini","Age","Mini","Mass","logL","logTe","logg","label","McoreTP",
         "C_O","period0","period1","pmode","Mloss","tau1m","X","Y","Xc","Xn",
         "Xo","Cexcess","Z","mbolmag","Gmag","G_BPmag","G_RPmag","B_Tmag",
         "V_Tmag","Jmag","Hmag","Ksmag"]
         
names = ["Zini","Age","Mini","Mass","logL","logTe","logg","label","McoreTP",
         "C_O","period0","period1","pmode","Mloss","tau1m","X","Y","Xc","Xn",
         "Xo","Cexcess","Z","mbolmag","Gmag","G_BPbrmag", "G_BPftmag",
         "G_RPmag","B_Tmag", "V_Tmag","Jmag","Hmag","Ksmag"]         
         

isochrones = pd.read_csv(isochrones_file, delim_whitespace=True, names=names, 
                         comment="#", dtype="float")
             
tois = pd.read_csv(toi_file, sep=",", header=0, comment="#", 
                   dtype={"GaiaDR2_ID":"string"})

# -----------------------------------------------------------------------------
# Compute Absolute Gmag
# -----------------------------------------------------------------------------
# Compute absolute magnitude
tois["Dist"] = 1000 / tois["Plx"]
tois["e_Dist"] = np.abs(-1000 * tois["Plx"]**-2 *tois["e_Plx"])
tois["abs_Gmag"] = tois["Gmag"] - 5*np.log10(tois["Dist"]/10)

tois["e_abs_Gmag"] = np.sqrt( (tois["Gmag"]*tois["e_Gmag"])**2 
                        + (-5*np.log10(tois["Dist"]/10))**2 
                        * np.abs(-5*tois["e_Dist"]/tois["Dist"]/np.log(10))**2)

ages = list(set(isochrones["Age"]))
ages.sort()

metallicities = list(set(isochrones["Zini"]))
metallicities.sort()

# Calculate colour
isochrones["BpRp"] = isochrones["G_BPbrmag"] - isochrones["G_RPmag"]


# -----------------------------------------------------------------------------
# Estimating masses
# -----------------------------------------------------------------------------
# ~5 Gyr age
age = ages[-5]
subset = isochrones.loc[isochrones["Age"]==age][isochrones["BpRp"]>=0.8][isochrones["Gmag"]>=4]
mass_est = interp1d(subset["BpRp"].values, subset["Mass"].values, kind="cubic") 

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
plt.close("all")
fig = plt.figure()
ax1 = fig.add_subplot(111)

# Plot Isochrones
for age in ages[14:]:
    ax1.plot(isochrones.loc[isochrones["Age"]==age]["BpRp"][:-3], 
             isochrones.loc[isochrones["Age"]==age]["Gmag"][:-3], "-", 
             label="%0.3f Gyr" % (age/10**9))
"""
for zzz in metallicities[::10]:
    ax1.plot(isochrones.loc[isochrones["Zini"]==zzz]["BpRp"][:-3], 
             isochrones.loc[isochrones["Zini"]==zzz]["Gmag"][:-3], "-", 
             label="Z = %0.8f" % zzz)#(zzz/10**9))
"""
# Plot each star
for star_i, star in tois.iterrows():
    if star["BpRp"] > 1.0:
        ax1.errorbar(star["BpRp"], star["abs_Gmag"], yerr=star["e_abs_Gmag"], 
                     fmt="r*")
        
        text = "%s (d=%0.1f pc)" % (star["GaiaDR2_ID"], star["Dist"])
        ax1.text(star["BpRp"], star["abs_Gmag"]-0.1, text, fontsize="5", 
                 color="blue", horizontalalignment='center')
                 
        # Estimate the mass
        try:
            mass = mass_est(star["BpRp"])
            print("%s --> %f M_sun" % (star["GaiaDR2_ID"], mass))
        except:
            pass

ax1.set_ylim(13, 1)
ax1.set_xlim(0.8, 3.5)
ax1.set_xlabel(r"$(B_P-R_P)$")
ax1.set_ylabel("Absolute Gmag")

plt.grid()
plt.legend(loc="best", prop={'size': 8})
plt.gcf().set_size_inches(16, 9)
plt.savefig("toi_isochrones.pdf")
