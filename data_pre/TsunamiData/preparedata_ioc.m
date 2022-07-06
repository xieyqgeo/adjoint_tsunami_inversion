clear;
close all;
lat0 = -57; lon0 = -25;
[StaInfo, nStations] = readStationsInfoCTL('iocStations.ctl');
% tw_removeTidesPolynomialFit('iocStations.ctl');   % remove tide using poly fit
dt = 10;          % sampling interval
tlen = 24000;     % total time, in seconds
t = 0:dt:tlen;          % desired data time vector
nt=length(t);
nsta=nStations;
data_obs=zeros(nsta,nt); 
% filter parameters
fl=1/5000; fh=1/200; 
sr=1/dt;
[BB,AA] = butter(2,[fl,fh]/(sr/2),'bandpass');
figure(1001); hold on;
for i=1:nsta
    fname=['./',StaInfo(i).name,'.dat'];
    tempdata=load(fname);
    tempdata_interp=interp1(tempdata(:,1),tempdata(:,2),t);
    tempdata1=filtfilt(BB,AA,tempdata_interp);
    data_obs(i,:)=tempdata1;
    plot(t/60,data_obs(i,:)/max(data_obs(i,:))+i,'b');
    text(0,i+0.5,StaInfo(i).name);
end
xlabel('time (min)'); title('origin tsunami recording (normalized amplitude)')
save data_obs_ori.mat data_obs StaInfo dt t

% plot the station
figure(1002); hold on;
worldmap('world');load coast; plotm(lat,long);
for i=1:nsta
    plotm(StaInfo(i).lat,StaInfo(i).lon,'^','markersize',10,'markeredgecolor','k')
end
plotm(lat0,lon0,'p','markersize',10,'markeredgecolor','r')


