# adjoint_tsunami_inversion
# tsunamiadjoint
A code package to the adjoint-state tsunami source inversion method
------
Author: Tong Zhou, Yuqing Xie

Acknowledgement: Chao An (The author of pcomcot, tsnuami data downloading codes) 

Dongdong Tian (The author of HinetPy if you use s-net data)
 
------
## 0. Theoretical background
For information see Zhou et al., (2019). An adjoint-state Full-waveform Tsunami Source Inversion Method and Its Application to the 2014 Chile-Iquique Tsunami Event. JGR: Solid Earth, 124(7), 6737-6750. doi:10.1029/2018JB016678

## 1. Prerequisities
MATLAB (recommend: R2016a or newer release)
Python3 (recommend: 3.6+ if you use s-net data)
Fortran compilier (gnu gfortran or intel ifort)

## 2. build package
### 2.1 build pcomcot, pcomcot_timereversal, pcomcot_fault
go to each folder and type 

`make nocdf`

Change compiler
### 2.2 build auxiliary fortran codes
go to `aux/` and type 

`make`

## 3. pre-process tsunami data
### 3.1 tsunami data download
go to `data_pre/TsunamiData`

use `iocDataDownload.py` and `DARTDataDownload.py`

change the download time span and station name. For DART informations, go to <u>https://nctr.pmel.noaa.gov/Dart/</u>. For ioc gauge informations, <u>http://www.ioc-sealevelmonitoring.org</u>

Downloaded file will be [StationName].txt and [ioc/DART]Stations.ctl (station information file). Originally, the sampling rate for ioc gauge is 60 seconds, and for DART there's dynamic sampling rate (Tsunami mode, 15 s).

The format of [StationName].txt  is:
T   h
…  …

  
 The format of Stations.ctl is:

Lon  lat name
  …  …   …


### 3.2 data processing
step1: go to `data_pre/TsunamiData`

step2: use `preparedata_ioc.m` [for ioc gauge data, first need to remove tide signal by polynomial fitting using 'tw_removeTidesPolynomialFit.m'.]
set sampling interval (dt) and time for recording (tlen). make sure to check the data (if there's NaNs)
output will be `data_obs_ori.mat` and [StationName].dat

then use `pick_first_wiggle.m` to pick up the first wiggle of tsunami waveform

If use synthetic data and use the whole waveforms without picking arrival time, the script preparedata_syn_withoutpicking_downsampling.m can be used. The input data data_obs_ori.mat saved a variable data_obs with the format:
data_obs=[ a1 a2 a3 a4…
B1 b2 b3 b4 …
…];
where a1, a2, a3 a3 … is the water level history of station a, b1, b2, b3, b4 are the water level history of station b
the output of step2 will be `data_obs_mod.mat`

step 3 use `create_timereversal_data.m` to create the data used in time reversal
the output will be `data_obs.mat` Stations.ctl and Station[XXXX].dat in `data_pre` folder


## 4. get bathymetry data
We suggest use GEBCO global grid (15 arc second resolution). For more information, see <u>https://www.gebco.net/data_and_products/gridded_bathymetry_data/</u>.
rename the downloaded file as `GEBCO_2020.nc`
go to `bathy_pre` and use `read_gebco.m` to cut the bathymetry data in the region covering all the stations and the source.
The three output files are _xcoordinate01.dat, _ycoordinate01.dat and Layer01.xyz.

_xcoordinate01.dat, _ycoordinate01.dat are row vectors，latitude is from south to north, longitude is from west to east  
Layer01.xyz: loop for latitude first (from south to north first, then from west to east).
An example:
126.004 25.0042 141
126.021 25.0042 124
...


## 5. running adjoint tsunami inversion

run `init_inversion.bash`. Change TsuData beforehead (if you use TsunamiData (for DART and Gauge) or S-net)

### 5.1 change parameters in the file `inversion_parameter.m`
simulation region: the simulation region for tsunami modelling. match the layer01.xyz (bathymetry file)

tsunami source box range: the region will be updated in the adjoint inversion

min and max water elevation: if inverted water elevation out of this range, it will be rejected

nsta and wt: number of stations and the weighting of each station. Change accordingly

other parameters read the comments in the file


### 5.2 run the main inversion program
`matlab -nodisplay -r driver_comcot`

use `screen` for offline running.


## 6 post process

download all the `inversion/` folder and use the `postprocess.m`.


## 7 process s-net data

For mac users, detailed instruction are in the `data_pre/S-net/s-net_mac_tutorial.pdf`

### 7.1 convert s-net data to sac and merge

Install HinetPy <u>https://seisman.github.io/HinetPy/installation.html</u>, compile win32tools and add the win32tools binaries in the environment variable

Put downloaded s-net data in `ori/`

use `extractsac_merge.py` to convert s-net data and merge. Output will be multiple sac files in `ori/`

### 7.2 convert sac to data file

use `read2ret.m`

### 7.3 convert data file to water elevation recordings useable in comcot

use `preparedata_snet.m`

#### 7.4 pick first wiggles

the same as section 3.2. use `pick_first_wiggle.m`

output data will be redirected to the upper folder (`data_pre/`)



## 8. running the pcomcot forward simulation
run `init_inversion.bash`. Change TsuData beforehead (if you use TsunamiData (for DART and Gauge) or S-net)


### 8.1 run pcomcot with fault rupture
goto `pcomcot_fault/`

modify `pcomcot.ctl`

modify `FaultParameters.ctl`. suggest to use centroid moment tensor solution (global CMT, <u>https://www.globalcmt.org</u>) and use scaling relation to estimate rupture length (L) and rupture width (W).

check `Stations.ctl`

run `mpirun -np X pcomcot`   X is desired cores you want to use

### 8.2 check output and tsunami station location calibration
Due to the inaccuracy of bathymetry, there might be some tide gauge stations on land, resulting in 0 in the output. use `check_waveform.m` to examine the data. use `COMCOT_readBinaryDataSnapshot.m` to check the wave snapshot (usually the last frame) to locate if all your stations are located in the wet grids (with non-zero wavefields). Move your stations on dry grids to the nearest wet grid and modify the Stations.ctl.

*Take notes on how far you moved the station. There might need a calibration (not implemented yet but will do in the future...)*

If modified `Stations.ctl`, replace it with the one in folder `data_pre/` and rerun `init_inversion.bash`
