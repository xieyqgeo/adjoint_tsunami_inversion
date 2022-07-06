% create_comcot_input
function create_comcot_input_reversal(datadir,src_grd,staweight)
% generate i_Timereversal, input: comcotdir, src_grd (in degree), station
% weighting
fid=fopen(['./i_TimeReversal'],'w');
fprintf(fid,'%s\n',datadir);   % data dir (contains station data for time reversal)
fprintf(fid,'%s\n','5');             % TR_NR
fprintf(fid,'%s\n','2');             % TR_PointSourceType
fprintf(fid,'%s\n','1');             % TR_AmpCorrection
fprintf(fid,'%s %s\n',num2str(src_grd(1)),num2str(src_grd(2)));      % TR_Epix, TR_Epiy
for i=1:length(staweight)-1
    fprintf(fid,'%s ',num2str(staweight(i)));
end
fprintf(fid,'%s',num2str(staweight(end)));
fclose(fid);
end
