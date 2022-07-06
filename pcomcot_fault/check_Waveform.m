% read and compare the selected taunami waves
clear;
close all;
load data_obs.mat
scale = 1;

nsta = length(StaInfo);
nt = t(end)+1;
dt = dt_compute;

data_syn = read_fwd(nt,dt,nsta,0,0,0);

obs_scale = max(max(data_obs));

figure(1901); hold on;
for i=1:nsta
    plot(t/60,data_obs_ori(i,:)/max(data_obs(i,:))+i*3,'k','linewidth',0.5);
    plot(t/60,data_obs(i,:)/max(data_obs(i,:))+i*3,'b','linewidth',1);
    plot(t(1:end)/60,data_syn(i,1:10:end)/max(data_syn(i,:))+i*3,'r','linewidth',1);
    %plot(t(1:end-1)/60,data_syn(i,:)/obs_scale*1e12+i,'r','linewidth',1);
    text(10,i*3+0.5,StaInfo(i).name);
end
xlabel('time (min)');