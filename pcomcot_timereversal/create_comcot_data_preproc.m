% create_comcot_input
function create_comcot_data_preproc(comcotdir,tmin,tmax,dt)
% generate para_for_conv, input: comcotdir, t_min, t_max, dt
% weighting
fid=fopen([comcotdir,'/Forward/para_for_conv'],'w');
fprintf(fid,'%s%s ',num2str(tmin),'d0');
fprintf(fid,'%s%s ',num2str(tmax),'d0');
fprintf(fid,'%s%s ',num2str(dt),'d0');
fclose(fid);
end


%mexico 0.0 17400.0 1.0
