close all;
clear all;

% using the whole waveform within the first 60 minutes

load data_obs_ori.mat;
[StaInfo, nStations] = readStationsInfoCTL('Stations.ctl');

dt=1;% new sampling rate, in second
total_time_min=60; % total time in minutes for inversion
dt_ori=0.5; % original sampling rate,in second

t=0:dt:total_time_min*60;
t=t(1:end-1);
data_obs2=data_obs;
for i=1:size(data_obs,1)
data_obs2(i,1:length(t))=interp1(0:dt_ori:size(data_obs,2)*dt_ori-dt_ori,data_obs(i,1:end),t);
end
data_obs=data_obs2(:,1:length(t));
t1=[zeros(size(data_obs,1),1) total_time_min*60+zeros(size(data_obs,1),1)];
t2=ones(size(data_obs,1),1)*total_time_min*60;
save data_obs_mod.mat data_obs StaInfo dt t t1 t2;
    

