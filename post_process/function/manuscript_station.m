% read station data
close all;
clear all;
nline=8;nrow=5;
max_t=120;
  ha=tight_subplot(nline,nrow,[0.04 0.04],[0.04 0.04],[0.04 0.04]);
k=0;
for i=1:38
     k=k+1;
    axes(ha(k));
    fn=['Station',num2str(i,'%04d'),'.dat']
station=COMCOT_readBinaryDataStation(fn);

plot(station(:,1)/60,station(:,2))
xlim([0 250])
end