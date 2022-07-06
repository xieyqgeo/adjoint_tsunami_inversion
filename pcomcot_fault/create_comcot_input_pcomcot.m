function create_comcot_input_pcomcot(init_cond,total_time,time_step,time_save)
% create pcomcot.ctl for time reversal
fid=fopen(['./pcomcot.ctl'],'w');
fprintf(fid,'%s\n','################################################################################');
fprintf(fid,'%s\n','#                                                                              #');             % TR_NR
fprintf(fid,'%s\n','#               Control file for PCOMCOT program (v F3.0 )                     #');
fprintf(fid,'%s\n','#                                                                              #');
fprintf(fid,'%s\n','################################################################################');
fprintf(fid,'%s\n','#===============================================:===============================');
fprintf(fid,'%s\n','# Parameters for Simulation                     :     Values                   |');
fprintf(fid,'%s\n','#===============================================:===============================');
fprintf(fid,'%s\n',' Coordinate System    (0:spherical, 1:cartesian):      0 ');
fprintf(fid,'%s\n',' Governing Equations  (0:linear,    1:nonlinear):      0 ');
fprintf(fid,'%s\n',' Specify Min WaterDepth Offshore         (meter):      0.00 ');
fprintf(fid,'%s %s\n',' Initial Condition             (0:file, 1:fault):    ',num2str(init_cond));
fprintf(fid,'%s\n',' Boundary Condition   (0:Open; 1:Wall; 2:Sponge):      1 ');
fprintf(fid,'%s %f\n',' Total run time                         (second):  ',total_time);
fprintf(fid,'%s %f\n',' Time step                              (second):   ',time_step);
fprintf(fid,'%s %f\n',' Time interval to Save Data             (second):  ',time_save);
fprintf(fid,'%s\n',' Save Flux                         (0:no, 1:yes):      0');
fprintf(fid,'%s\n',' Compute Green s Functions         (0:no, 1:yes):      0');
fprintf(fid,'%s\n',' Minimum grids on each computing node           :   2000');
fprintf(fid,'%s\n',' Feedback to parent layer          (0:no, 1:yes):      0');
fprintf(fid,'\n');
fprintf(fid,'%s\n','#===============================================:===============================');
fprintf(fid,'%s\n','# System Parameters (DO NOT CHANGE)             :     Values                   |');
fprintf(fid,'%s\n','#===============================================:===============================');
fprintf(fid,'%s\n',' Initial Water Elevation File Name              :  InitialWaterElevation.nf/xyz');
fprintf(fid,'%s\n',' File Name of Bathymetry Data                   :  layerXX.nf/xyz');
fprintf(fid,'%s\n',' File Name of Stations                          :  Stations.txt');
fclose(fid);
end
