#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import numpy.matlib
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import datetime
import netCDF4 as nc
from scipy.stats import linregress
import copy


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def ubrmse(predictions, targets):
    return np.sqrt((((predictions - np.nanmean(predictions)) -
                    (targets - np.nanmean(targets)))**2).mean())


if len(sys.argv) > 1:
    cname = str(sys.argv[1])
else:
    cname = "BE-Bra"    

layers = np.array([1, 3, 6])
layer_names = ["5cm", "20cm", "50cm"]

startyear = 2009
endyear = 2018
nyears = endyear - startyear + 1

startdate = datetime.datetime(startyear, 1, 1)
enddate = datetime.datetime(endyear, 12, 31)
timedelta = dates.num2timedelta(1.0 / 24.0) # 1 hour = 1 / 24 of a day
xdate = dates.drange(startdate, enddate, timedelta)
yl = dates.YearLocator()
ml = dates.MonthLocator()
yearsFmt = dates.DateFormatter("%Y")
monthsFmt = dates.DateFormatter('%b')


OBS_mean = np.zeros((nyears, 365, 3))
OL_mean_swc = np.zeros((nyears, 365, 3))
DAS_mean_swc = np.zeros((nyears, 365, 3))
DASP_mean_swc = np.zeros((nyears, 365, 3))
DASHP_mean_swc = np.zeros((nyears, 365, 3))

OBS_mean_et = np.zeros((nyears, 365))
OL_mean_et = np.zeros((nyears, 365))
DAS_mean_et = np.zeros((nyears, 365))
DASP_mean_et = np.zeros((nyears, 365))
DASHP_mean_et = np.zeros((nyears, 365))

params_mean = np.zeros((nyears, 365, 3, 4))
params_std = np.zeros((nyears, 365, 3, 4))

day_i = 0

for yyyy in range(startyear, endyear+1):
    y = yyyy - startyear
    OBS_path = "/p/project/cjicg41/jicg4177/Observations/" + cname + "_swc.nc"
    OL_path = cname + "-Ens/Ensemble_Collection_OL."+str(yyyy)+".nc"
    DAS_path = cname + "-Ens/Ensemble_Collection_DAS."+str(yyyy)+".nc"
    DASP_path = cname + "-Ens/Ensemble_Collection_DASP."+str(yyyy)+".nc"
    DASHP_path = cname + "-Ens/Ensemble_Collection_DASHP."+str(yyyy)+".nc"

    OBS_src = nc.Dataset(OBS_path, "r")
    OL_src = nc.Dataset(OL_path, "r")
    DAS_src = nc.Dataset(DAS_path, "r")
    DASP_src = nc.Dataset(DASP_path, "r")
    DASHP_src = nc.Dataset(DASHP_path, "r")

    OBS_mean_swc = np.reshape(OBS_src["obs_swc"][:], (11, 365, 3))
    OL_mean = OL_src.variables["H2OSOI_mean"][:]  * 100.0
    DAS_mean = DAS_src.variables["H2OSOI_mean"][:]  * 100.0
    DASP_mean = DASP_src.variables["H2OSOI_mean"][:] * 100.0
    DASHP_mean = DASHP_src.variables["H2OSOI_mean"][:] * 100.0

#    OL_mean = OL_src.variables["SOILLIQ_mean"][:]  / 1000.0 / 0.04 * 100.0
#    DAS_mean = DAS_src.variables["SOILLIQ_mean"][:]  / 1000.0 / 0.04 * 100.0

    OL_mean_swc[y, :, 0] = OL_mean[:, layers[0], 0]
    OL_mean_swc[y, :, 1] = (OL_mean[:, layers[1], 0] + OL_mean[:, layers[1]+1, 0]) / 2.0
    OL_mean_swc[y, :, 2] = OL_mean[:, layers[2], 0]

    DAS_mean_swc[y, :, 0] = DAS_mean[:, layers[0], 0]
    DAS_mean_swc[y, :, 1] = (DAS_mean[:, layers[1], 0] + DAS_mean[:, layers[1]+1, 0]) / 2.0
    DAS_mean_swc[y, :, 2] = DAS_mean[:, layers[2], 0]

    DASP_mean_swc[y, :, 0] = DASP_mean[:, layers[0], 0]
    DASP_mean_swc[y, :, 1] = (DASP_mean[:, layers[1], 0] + DASP_mean[:, layers[1]+1, 0]) / 2.0
    DASP_mean_swc[y, :, 2] = DASP_mean[:, layers[2], 0]

    DASHP_mean_swc[y, :, 0] = DASHP_mean[:, layers[0], 0]
    DASHP_mean_swc[y, :, 1] = (DASHP_mean[:, layers[1], 0] + DASHP_mean[:, layers[1]+1, 0]) / 2.0
    DASHP_mean_swc[y, :, 2] = DASHP_mean[:, layers[2], 0]

    try:
        OBS_mean_et = np.reshape(OBS_src["obs_et"][:], (11, 365))
    except IndexError:
        print("ET not available in " + str(OBS_path) + " set to 0.")
        OBS_mean_et = np.zeros((11, 365))
    OL_et = OL_src.variables["QFLX_EVAP_TOT_mean"][:] * 86400.0
    DAS_et = DAS_src.variables["QFLX_EVAP_TOT_mean"][:] * 86400.0
    DASP_et = DASP_src.variables["QFLX_EVAP_TOT_mean"][:] * 86400.0
    DASHP_et = DASHP_src.variables["QFLX_EVAP_TOT_mean"][:] * 86400.0


    OL_mean_et[y, :] = OL_et[:,0]
    DAS_mean_et[y, :] = DAS_et[:,0]
    DASP_mean_et[y,:] = DASP_et[:,0]
    DASHP_mean_et[y,:] = DASHP_et[:,0]
    for l in range(3):
        params_mean[y, :, l, 0] = DASHP_src.variables["DA_bsw_mean"][:, layers[l], 0]
        params_mean[y, :, l, 1] = DASHP_src.variables["DA_sucsat_mean"][:, layers[l], 0]
        params_mean[y, :, l, 2] = DASHP_src.variables["DA_hksat_min_mean"][:, layers[l], 0]
        params_mean[y, :, l, 3] = DASHP_src.variables["DA_watsat_mean"][:, layers[0], 0]
        params_std[y, :, l, 0] = DASHP_src.variables["DA_bsw_std"][:, layers[l], 0]
        params_std[y, :, l, 1] = DASHP_src.variables["DA_sucsat_std"][:, layers[l], 0]
        params_std[y, :, l, 2] = DASHP_src.variables["DA_hksat_min_std"][:, layers[l], 0]
        params_std[y, :, l, 3] = DASHP_src.variables["DA_watsat_std"][:, layers[0], 0]



# PLOT 1 : time series
fig = plt.figure(figsize=(16,9))
plt.rcParams.update({'font.size': 20})
axlist = [fig.add_subplot(311), fig.add_subplot(312), fig.add_subplot(313)]
date_days = xdate[0:-1:24]

for l in range(3):
    axlist[l].plot(date_days[:-1], OBS_mean_swc[:-1, :, l].flatten(), color="darkred", label="OBS")
    axlist[l].plot(date_days[:-1], OL_mean_swc[:, :, l].flatten(), color="green", label="OL")
    axlist[l].plot(date_days[:-1], DAS_mean_swc[:, :, l].flatten(), color="blue", label="DA_s")
    axlist[l].plot(date_days[:-1], DASP_mean_swc[:, :, l].flatten(), color="teal", label="DA s+p+o")
    axlist[l].plot(date_days[:-1], DASHP_mean_swc[:, :, l].flatten(), color="purple", label="DA s+hp")

    axlist[l].set_ylim(0.0, 70.0)
    axlist[l].set_xlim(date_days[0] - 30, date_days[-1] + 30)
    axlist[l].yaxis.set_label_coords(-0.08, 0.4)
    axlist[l].yaxis.set_major_locator(ticker.MultipleLocator(15.0))
    axlist[l].xaxis.set_major_locator(yl)
    axlist[l].xaxis.set_major_formatter(yearsFmt)
    axlist[l].xaxis.set_minor_locator(ml)
    if l == 2:
        axlist[l].xaxis.set_major_formatter(yearsFmt)
    else:
        axlist[l].xaxis.set_ticklabels([])

axlist[0].set_title("Depth 5cm")
axlist[1].set_title("Depth 20cm")
axlist[2].set_title("Depth 50cm")

axlist[0].set_ylabel("SWC [%]", rotation="horizontal", position=(0.005,1.05))
axlist[1].set_ylabel("SWC [%]", rotation="horizontal", position=(0.005,1.05))
axlist[2].set_ylabel("SWC [%]", rotation="horizontal", position=(0.005,1.05))

handles, labels = axlist[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', fontsize=8, ncol=5)
fig.tight_layout()
plt.subplots_adjust(bottom=0.11)
outfile = "Plot_Timeseries_" + str(cname)
plt.savefig(outfile + ".png", format="png", dpi=600)


#PLOT 2: correlations
for l in range(3):
    fig = plt.figure(figsize=(16,16))
    plt.rcParams.update({'font.size': 20})
    axlist = [fig.add_subplot(221), fig.add_subplot(222), fig.add_subplot(223), fig.add_subplot(224)]
    date_days = xdate[0:-1:24]
    lw=1.0
    
    axlist[0].scatter(OBS_mean_swc[:-1, :, l].flatten(), OL_mean_swc[:, :, l].flatten(), marker="x", s=30, color="green", label="OL", linewidth=lw, alpha=0.8)
    
    axlist[1].scatter(OBS_mean_swc[:-1, :, l].flatten(), DAS_mean_swc[:, :, l].flatten(), marker="o", s=30, color="blue", label="DA s", linewidth=lw, alpha=0.8)
    
    axlist[2].scatter(OBS_mean_swc[:-1, :, l].flatten(), DASP_mean_swc[:, :, l].flatten(), marker="+", s=30, color="teal", label="DA s+p+o", linewidth=lw, alpha=0.8)
    
    axlist[3].scatter(OBS_mean_swc[:-1, :, l].flatten(), DASHP_mean_swc[:, :, l].flatten(), marker="v", s=30, color="purple", label="DA s+hp", linewidth=lw, alpha=0.8)
    
    
    axlist[0].plot(np.arange(5.0, 70.0, 5.0),np.arange(5.0, 70.0, 5.0), linewidth=1.5, color="black")
    axlist[1].plot(np.arange(5.0, 70.0, 5.0),np.arange(5.0, 70.0, 5.0), linewidth=1.5, color="black")
    axlist[2].plot(np.arange(5.0, 70.0, 5.0),np.arange(5.0, 70.0, 5.0), linewidth=1.5, color="black")
    axlist[3].plot(np.arange(5.0, 70.0, 5.0),np.arange(5.0, 70.0, 5.0), linewidth=1.5, color="black")
    
    axlist[0].set_ylabel("OL - SWC [%]", rotation="horizontal")
    axlist[1].set_ylabel("DA s - SWC [%]", rotation="horizontal")
    axlist[2].set_ylabel("DA s+p+o - SWC [%]", rotation="horizontal")
    axlist[3].set_ylabel("DA s+hp - SWC [%]", rotation="horizontal")
    
    axlist[0].set_xlabel("OBS", rotation="horizontal")
    axlist[1].set_xlabel("OBS", rotation="horizontal")
    axlist[2].set_xlabel("OBS", rotation="horizontal")
    axlist[3].set_xlabel("OBS", rotation="horizontal")
    
    axlist[0].get_yaxis().set_label_coords(0.01,1.01)
    axlist[1].get_yaxis().set_label_coords(0.01,1.01)
    axlist[2].get_yaxis().set_label_coords(0.01,1.01)
    axlist[3].get_yaxis().set_label_coords(0.01,1.01)
    
    axlist[0].xaxis.set_major_locator(ticker.MultipleLocator(10.0))
    axlist[0].yaxis.set_major_locator(ticker.MultipleLocator(10.0))
    axlist[1].xaxis.set_major_locator(ticker.MultipleLocator(10.0))
    axlist[1].yaxis.set_major_locator(ticker.MultipleLocator(10.0))
    axlist[2].xaxis.set_major_locator(ticker.MultipleLocator(10.0))
    axlist[2].yaxis.set_major_locator(ticker.MultipleLocator(10.0))
    axlist[3].xaxis.set_major_locator(ticker.MultipleLocator(10.0))
    axlist[3].yaxis.set_major_locator(ticker.MultipleLocator(10.0))
    
    textx = 42
    texty1 = 9
    texty2 = 14
    texty3 = 19
    texty4 = 24
    
    OBS_reduced = OBS_mean_swc[:-1, :, l].flatten()[~np.isnan(OBS_mean_swc[:-1, :, l].flatten())]
    OL_reduced = OL_mean_swc[:, :, l].flatten()[~np.isnan(OBS_mean_swc[:-1, :, l].flatten())]
    DAS_reduced = DAS_mean_swc[:, :, l].flatten()[~np.isnan(OBS_mean_swc[:-1, :, l].flatten())]
    DASP_reduced = DASP_mean_swc[:, :, l].flatten()[~np.isnan(OBS_mean_swc[:-1, :, l].flatten())]
    DASHP_reduced = DASHP_mean_swc[:, :, l].flatten()[~np.isnan(OBS_mean_swc[:-1, :, l].flatten())]
    
    axlist[0].text(textx, texty1, r"$R^2$ = {0:.2f}".format(linregress(OBS_reduced, OL_reduced).rvalue),
              fontsize=20)
    axlist[1].text(textx, texty1, r"$R^2$ = {0:.2f}".format(linregress(OBS_reduced, DAS_reduced).rvalue),
              fontsize=20)
    axlist[2].text(textx, texty1, r"$R^2$ = {0:.2f}".format(linregress(OBS_reduced, DASP_reduced).rvalue),
              fontsize=20)
    axlist[3].text(textx, texty1, r"$R^2$ = {0:.2f}".format(linregress(OBS_reduced, DASHP_reduced).rvalue),
              fontsize=20)
    
    axlist[0].text(textx, texty4, "RMSE = {0:.2f}".format(rmse(OL_reduced, OBS_reduced)),
              fontsize=20)
    axlist[1].text(textx, texty4, "RMSE = {0:.2f}".format(rmse(DAS_reduced, OBS_reduced)),
              fontsize=20)
    axlist[2].text(textx, texty4, "RMSE = {0:.2f}".format(rmse(DASP_reduced, OBS_reduced)),
              fontsize=20)
    axlist[3].text(textx, texty4, "RMSE = {0:.2f}".format(rmse(DASHP_reduced, OBS_reduced)),
              fontsize=20)
    
    
    axlist[0].text(textx, texty2, "MBE = {0:.2f}".format(np.mean((-OL_reduced + OBS_reduced))),
              fontsize=20)
    axlist[1].text(textx, texty2, "MBE = {0:.2f}".format(np.mean((-OL_reduced + DAS_reduced))),
              fontsize=20)
    axlist[2].text(textx, texty2, "MBE = {0:.2f}".format(np.mean((-OL_reduced + DASP_reduced))),
              fontsize=20)
    axlist[3].text(textx, texty2, "MBE = {0:.2f}".format(np.mean((-OL_reduced + DASHP_reduced))),
              fontsize=20)
    
    
    axlist[0].text(textx, texty3, "ubRMSE = {0:.2f}".format(ubrmse(OL_reduced, OBS_reduced)),
              fontsize=20)
    axlist[1].text(textx, texty3, "ubRMSE = {0:.2f}".format(ubrmse(DAS_reduced, OBS_reduced)),
              fontsize=20)
    axlist[2].text(textx, texty3, "ubRMSE = {0:.2f}".format(ubrmse(DASP_reduced, OBS_reduced)),
              fontsize=20)
    axlist[3].text(textx, texty3, "ubRMSE = {0:.2f}".format(ubrmse(DASHP_reduced, OBS_reduced)),
              fontsize=20)
    
    axlist[0].set_aspect("equal")
    axlist[1].set_aspect("equal")
    axlist[2].set_aspect("equal")
    axlist[3].set_aspect("equal")
    
    fig.tight_layout()
    outfile = "Plot_CLM_Scatter_" + str(cname) + "_" + layer_names[l]
    plt.savefig(outfile + ".png", format="png", dpi=300)


# PLOT 3 : time series ET
fig = plt.figure(figsize=(16,9))
plt.rcParams.update({'font.size': 20})
axlist = [fig.add_subplot(111)]
date_days = xdate[0:-1:24]
l = 0
axlist[l].plot(date_days[:-1], OBS_mean_et[:-1, :].flatten(), color="darkred", label="OBS")
axlist[l].plot(date_days[:-1], OL_mean_et[:, :].flatten(), color="green", label="OL")
axlist[l].plot(date_days[:-1], DAS_mean_et[:, :].flatten(), color="blue", label="DA_s")
axlist[l].plot(date_days[:-1], DASP_mean_et[:, :].flatten(), color="teal", label="DA s+p+o")
axlist[l].plot(date_days[:-1], DASHP_mean_et[:, :].flatten(), color="purple", label="DA s+hp")

#axlist[l].set_ylim(0.0, 70.0)
axlist[l].set_xlim(date_days[0] - 30, date_days[-1] + 30)
axlist[l].yaxis.set_label_coords(-0.08, 0.4)
axlist[l].yaxis.set_major_locator(ticker.MultipleLocator(1.0))
axlist[l].xaxis.set_major_locator(yl)
axlist[l].xaxis.set_major_formatter(yearsFmt)
axlist[l].xaxis.set_minor_locator(ml)
axlist[l].xaxis.set_major_formatter(yearsFmt)

axlist[0].set_ylabel("ET [mm/d]", rotation="horizontal", position=(0.005,1.05))

handles, labels = axlist[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', fontsize=8, ncol=5)
fig.tight_layout()
plt.subplots_adjust(bottom=0.11)
outfile = "Plot_Timeseries_ET_" + str(cname)
plt.savefig(outfile + ".png", format="png", dpi=600)


# PLOT 4 : scatter ET

fig = plt.figure(figsize=(16,16))
plt.rcParams.update({'font.size': 20})
axlist = [fig.add_subplot(221), fig.add_subplot(222), fig.add_subplot(223), fig.add_subplot(224)]
date_days = xdate[0:-1:24]
lw=1.0

axlist[0].scatter(OBS_mean_et[:-1, :].flatten(), OL_mean_et[:, :].flatten(), marker="x", s=30, color="green", label="OL", linewidth=lw, alpha=0.8)

axlist[1].scatter(OBS_mean_et[:-1, :].flatten(), DAS_mean_et[:, :].flatten(), marker="o", s=30, color="blue", label="DA s", linewidth=lw, alpha=0.8)

axlist[2].scatter(OBS_mean_et[:-1, :].flatten(), DASP_mean_et[:, :].flatten(), marker="+", s=30, color="teal", label="DA s+p+o", linewidth=lw, alpha=0.8)

axlist[3].scatter(OBS_mean_et[:-1, :].flatten(), DASHP_mean_et[:, :].flatten(), marker="v", s=30, color="purple", label="DA s+hp", linewidth=lw, alpha=0.8)

axlist[0].plot(np.arange(-1.0, 11.0, 1.0),np.arange(-1.0, 11.0, 1.0), linewidth=1.5, color="black")
axlist[1].plot(np.arange(-1.0, 11.0, 1.0),np.arange(-1.0, 11.0, 1.0), linewidth=1.5, color="black")
axlist[2].plot(np.arange(-1.0, 11.0, 1.0),np.arange(-1.0, 11.0, 1.0), linewidth=1.5, color="black")
axlist[3].plot(np.arange(-1.0, 11.0, 1.0),np.arange(-1.0, 11.0, 1.0), linewidth=1.5, color="black")

axlist[0].set_ylabel("OL - ET [mm/d]", rotation="horizontal")
axlist[1].set_ylabel("DA s - ET [mm/d]", rotation="horizontal")
axlist[2].set_ylabel("DA s+p+o - ET [mm/d]", rotation="horizontal")
axlist[3].set_ylabel("DA s+hp - ET [mm/d]", rotation="horizontal")

axlist[0].set_xlabel("OBS", rotation="horizontal")
axlist[1].set_xlabel("OBS", rotation="horizontal")
axlist[2].set_xlabel("OBS", rotation="horizontal")
axlist[3].set_xlabel("OBS", rotation="horizontal")

axlist[0].get_yaxis().set_label_coords(0.01,1.01)
axlist[1].get_yaxis().set_label_coords(0.01,1.01)
axlist[2].get_yaxis().set_label_coords(0.01,1.01)
axlist[3].get_yaxis().set_label_coords(0.01,1.01)

axlist[0].xaxis.set_major_locator(ticker.MultipleLocator(1.0))
axlist[0].yaxis.set_major_locator(ticker.MultipleLocator(1.0))
axlist[1].xaxis.set_major_locator(ticker.MultipleLocator(1.0))
axlist[1].yaxis.set_major_locator(ticker.MultipleLocator(1.0))
axlist[2].xaxis.set_major_locator(ticker.MultipleLocator(1.0))
axlist[2].yaxis.set_major_locator(ticker.MultipleLocator(1.0))
axlist[3].xaxis.set_major_locator(ticker.MultipleLocator(1.0))
axlist[3].yaxis.set_major_locator(ticker.MultipleLocator(1.0))

axlist[0].set_ylim(-1.5, 11.5)
axlist[0].set_xlim(-1.5, 11.5)
axlist[1].set_ylim(-1.5, 11.5)
axlist[1].set_xlim(-1.5, 11.5)
axlist[2].set_ylim(-1.5, 11.5)
axlist[2].set_xlim(-1.5, 11.5)
axlist[3].set_ylim(-1.5, 11.5)
axlist[3].set_xlim(-1.5, 11.5)

textx = 6.7
texty1 = 0.1
texty2 = 1.1
texty3 = 2.1
texty4 = 3.1


OBS_reduced = OBS_mean_et[:-1, :].flatten()[~np.isnan(OBS_mean_et[:-1, :].flatten())]
OL_reduced = OL_mean_et[:, :].flatten()[~np.isnan(OBS_mean_et[:-1, :].flatten())]
DAS_reduced = DAS_mean_et[:, :].flatten()[~np.isnan(OBS_mean_et[:-1, :].flatten())]
DASP_reduced = DASP_mean_et[:, :].flatten()[~np.isnan(OBS_mean_et[:-1, :].flatten())]
DASHP_reduced = DASHP_mean_et[:, :].flatten()[~np.isnan(OBS_mean_et[:-1, :].flatten())]

axlist[0].text(textx, texty1, r"$R^2$ = {0:.2f}".format(linregress(OBS_reduced, OL_reduced).rvalue),
          fontsize=20)
axlist[1].text(textx, texty1, r"$R^2$ = {0:.2f}".format(linregress(OBS_reduced, DAS_reduced).rvalue),
          fontsize=20)
axlist[2].text(textx, texty1, r"$R^2$ = {0:.2f}".format(linregress(OBS_reduced, DASP_reduced).rvalue),
          fontsize=20)
axlist[3].text(textx, texty1, r"$R^2$ = {0:.2f}".format(linregress(OBS_reduced, DASHP_reduced).rvalue),
          fontsize=20)

axlist[0].text(textx, texty4, "RMSE = {0:.2f}".format(rmse(OL_reduced, OBS_reduced)),
          fontsize=20)
axlist[1].text(textx, texty4, "RMSE = {0:.2f}".format(rmse(DAS_reduced, OBS_reduced)),
          fontsize=20)
axlist[2].text(textx, texty4, "RMSE = {0:.2f}".format(rmse(DASP_reduced, OBS_reduced)),
          fontsize=20)
axlist[3].text(textx, texty4, "RMSE = {0:.2f}".format(rmse(DASHP_reduced, OBS_reduced)),
          fontsize=20)

axlist[0].text(textx, texty2, "MBE = {0:.2f}".format(np.mean((-OL_reduced + OBS_reduced))),
          fontsize=20)
axlist[1].text(textx, texty2, "MBE = {0:.2f}".format(np.mean((-OL_reduced + DAS_reduced))),
          fontsize=20)
axlist[2].text(textx, texty2, "MBE = {0:.2f}".format(np.mean((-OL_reduced + DASP_reduced))),
          fontsize=20)
axlist[3].text(textx, texty2, "MBE = {0:.2f}".format(np.mean((-OL_reduced + DASHP_reduced))),
          fontsize=20)

axlist[0].text(textx, texty3, "ubRMSE = {0:.2f}".format(ubrmse(OL_reduced, OBS_reduced)),
          fontsize=20)
axlist[1].text(textx, texty3, "ubRMSE = {0:.2f}".format(ubrmse(DAS_reduced, OBS_reduced)),
          fontsize=20)
axlist[2].text(textx, texty3, "ubRMSE = {0:.2f}".format(ubrmse(DASP_reduced, OBS_reduced)),
          fontsize=20)
axlist[3].text(textx, texty3, "ubRMSE = {0:.2f}".format(ubrmse(DASHP_reduced, OBS_reduced)),
          fontsize=20)

axlist[0].set_aspect("equal")
axlist[1].set_aspect("equal")
axlist[2].set_aspect("equal")
fig.tight_layout()
outfile = "Plot_CLM_ET_Scatter_" + str(cname)
plt.savefig(outfile + ".png", format="png", dpi=300)

