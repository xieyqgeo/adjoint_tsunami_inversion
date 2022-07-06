################################################################################
#                                                                              #
#               Control file for PCOMCOT program (v F3.0 )                     #
#                                                                              #
################################################################################
#===============================================:===============================
# Parameters for Simulation                     :     Values                   |
#===============================================:===============================
 Coordinate System    (0:spherical, 1:cartesian):      0 
 Governing Equations  (0:linear,    1:nonlinear):      0 
 Specify Min WaterDepth Offshore         (meter):      0.00 
 Initial Condition             (0:file, 1:fault):    1
 Boundary Condition   (0:Open; 1:Wall; 2:Sponge):      1 
 Total run time                         (second):   3600.000000
 Time step                              (second):    1.000000
 Time start to record snapshot          (second):    3400.000000
 Time interval to Save Data             (second):   200.000000
 Save Flux                         (0:no, 1:yes):      0
 Compute Green s Functions         (0:no, 1:yes):      0
 Minimum grids on each computing node           :   2000
 Feedback to parent layer          (0:no, 1:yes):      0

#===============================================:===============================
# System Parameters (DO NOT CHANGE)             :     Values                   |
#===============================================:===============================
 Initial Water Elevation File Name              :  InitialWaterElevation.nf/xyz
 File Name of Bathymetry Data                   :  layerXX.nf/xyz
 File Name of Stations                          :  Stations.txt
