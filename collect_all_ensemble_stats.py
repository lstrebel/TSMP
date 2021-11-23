#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import netCDF4 as nc
import datetime
import sys

if len(sys.argv) == 4:
    year = int(sys.argv[1])
    experiment = str(sys.argv[2])
    num_ensemble = int(sys.argv[3])
else:
    year = 2010
    experiment = "DA"
    num_ensemble = 96

if experiment == "OL":
    ensemble_description = "OpenLoop."
else:
    ensemble_description = "Data Assimilation SoilNet daily SWC."

ensemble_description += " Using perturbed soil properties (sand&clay % [uniform] and forcings (precipitation, short- and longwave radiation, 2m temperature [correlated according to ]."

if experiment == "OL":
    dst = nc.Dataset("Ensemble_Collection_OL." + str(year) + ".nc", "w", format="NETCDF4")
else:
    dst = nc.Dataset("Ensemble_Collection_DA." + str(year) + ".nc", "w", format="NETCDF4")
if experiment == "OL":
    cname = "WU_OL"
else:
    cname = "WU_DA"
base_src = nc.Dataset(str(cname)+".clm2_0000.h0."+str(year)+"-01-01-00000.nc", "r")

# copy attributes
for name in base_src.ncattrs():
    dst.setncattr("original_attribute_" + name, base_src.getncattr(name))
# copy / reduce dimensions
for name, dimension in base_src.dimensions.items():
        dst.createDimension( name, len(dimension))
dst.createDimension("num_ens", num_ensemble)
dst.setncattr("Collected_Ensemble_From", "CLM+PDAF (TSMP)")

ensemble_dimensions = [('time', 'lndgrid'),
                       ('time', 'levsoi',  'lndgrid'),
                       ('time', 'levgrnd', 'lndgrid')]

# Copy all variables not to be collected into ensemble statistics
for name, var in base_src.variables.items():
    if not var.dimensions in ensemble_dimensions:
        dst.createVariable(name, var.datatype, var.dimensions)
        dst[name].setncatts(base_src[name].__dict__)
        dst[name][:] = base_src[name][:]
#print("----- copied all non-ensemble variables -------")

# Collect ensemble variables into new file
for name, var in base_src.variables.items():
    if not var.dimensions in ensemble_dimensions:
        continue
    #print(name)
    # Get size of dimensions for the current variable
    dim_vals = tuple([base_src.dimensions[base_src.variables[name].dimensions[i]].size for i in range(len(base_src.variables[name].dimensions))])
    dim_names = tuple([base_src.dimensions[base_src.variables[name].dimensions[i]].name for i in range(len(base_src.variables[name].dimensions))])

    min_val = np.inf * np.ones(dim_vals)
    min_ind = -1 * np.ones(dim_vals)
    max_val = -1.0 * np.ones(dim_vals)
    max_ind = -1 * np.ones(dim_vals)
    mean_val = -1.0 * np.ones(dim_vals)
    new_mean = np.zeros(dim_vals)
    prevariance = np.zeros(dim_vals)
    std_val = np.zeros(dim_vals)

    # Iterate through ensembles
    for ens in range(0, num_ensemble):
        src = nc.Dataset(str(cname)+".clm2_" + str(ens).zfill(4) + ".h0."+str(year)+"-01-01-00000.nc", "r").variables[name][:]

        if ens == 0:
            min_val[:] = src.data[:]
            max_val[:] = src.data[:]
            mean_val[:] = src.data[:]
            min_ind[:] = 0
            max_ind[:] = 0
        else:
            min_val[:] = np.minimum(min_val[:], src.data[:])
            min_ind[np.where(min_val[:] == src.data[:])] = ens

            max_val[:] = np.maximum(max_val[:], src.data[:])
            max_ind[np.where(max_val[:] == src.data[:])] = ens

            new_mean[:] = np.add(mean_val[:], np.subtract(src.data[:], mean_val[:]) / (ens + 1))
            prevariance[:] = np.add(prevariance[:], np.multiply(np.subtract(src.data[:], mean_val[:]),np.subtract(src.data[:], new_mean[:])))
            mean_val[:] = new_mean[:]

    std_val[:] = np.sqrt(prevariance / (num_ensemble - 1))

    ensemble_names = ["min", "min_ind", "max", "max_ind", "mean", "std"]
    ensemble_var = [min_val, min_ind, max_val, max_ind, mean_val, std_val]
    for e in range(len(ensemble_names)):
        typ = np.float64 if ensemble_names[e][-3:] != "ind" else np.int32
        fval = np.nan if ensemble_names[e][-3:] != "ind" else -1
        e_var = dst.createVariable(name + "_" + ensemble_names[e], typ, dim_names, fill_value=fval)
        e_var[:] = ensemble_var[e][:]
        e_var.units = base_src.variables[name].units if ensemble_names[e][-3:] != "ind" else "ensemble index"
        e_var.long_name = base_src.variables[name].long_name + " " + ensemble_names[e]

dst.description = "Collected ensemble statistics from CLM3.5+PDAF (terrsysmp) setup of 2 CORDEX gridcells corresponding to the study site Wuestebach for the year " + str(year) + "." + ensemble_description
dst.history = "Created " + datetime.datetime.today().strftime("%d.%m.%y")
dst.close()
#print("--------- collected all ensemble variables -----------")
