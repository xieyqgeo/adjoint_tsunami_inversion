% convert_station_file_to_comcot
load data_obs_mod.mat
dt_interp = dt;
tmax = max(t)+dt;
nsta = length(StaInfo);

% generate cut para
fid = fopen('cut_sta_para.ctl','w');
fprintf(fid,'%s %s %s','0.0d0',[num2str(tmax),'.0d0'],[num2str(dt_interp),'.0d0']);
fclose(fid);

% generate stations.ctl
fid = fopen('Stations.ctl','w');

for i=1:nsta
    fprintf(fid,'%s %s %s\n',num2str(StaInfo(i).lon),num2str(StaInfo(i).lat),StaInfo(i).name);
    % generate dat2 files
    fid2 = fopen([StaInfo(i).name,'.dat2'],'w');
    for j=1:length(t)
        fprintf(fid2,'%s %s\n',num2str(t(j)),num2str(data_obs(i,j)));
    end
    fclose(fid2);
end
fclose(fid);    

% generate data_obs.mat
data_obs_cut = data_obs;
StaInfo_cut = StaInfo;
dt_compute = dt_interp;
load data_obs_ori.mat
for i=1:length(StaInfo_cut)
    idx = find(strcmp({StaInfo.name},StaInfo_cut(i).name));
    data_obs_ori(i,:) = data_obs(idx,:);
end
data_obs = data_obs_cut;
StaInfo = StaInfo_cut;
timewindow_start = t1;
timewindow_end   = t2;
save data_obs.mat data_obs data_obs_ori StaInfo t dt dt_compute timewindow_start timewindow_end

system(['cp *.dat2 ../']);
system(['cp Stations.ctl ../']);
system(['cp cut_sta_para.ctl ../']);
cd ..
system('./convertTsunamiDataTRInput');
system('cp ./TsunamiData/data_obs.mat ./');
