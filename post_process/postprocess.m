% post_process_ajoint
clear all;
close all;
%% change parameter for plotting
iter_plot=10; % plot the result of which iteration
nline=6;nrow=6; % lines and rows of the subplot for waveforms comparison of each station 
gap_contour=2; % interval between two contour lines, in meter
caxis_h=[-8 8]; % water eleavtion range for plotting the inverted source
result_folder='./../inversion'

%%
addpath('./function');
load([result_folder '/tsunami_inversion_result.mat']) % After the computation is completed, all the variables and parameters are saved to this file in the folder "inversion" 
xrange_plot=[srclonmn srclonmx]; % plot range of longitude
yrange_plot=[srclatmn srclatmx]; % plot range of latitude
dt=time_step;
t=0:dt:total_time;
sta_start=timewindow_start;
sta_end=timewindow_end;
stalon=src_grd(1);
stalat=src_grd(2);


%% plot station recordings
load([result_folder '/data_pred_' num2str(iter_plot) '.mat']);
tsu_pred3=tsu_pred;

n=0;
stalon=[];
stalat=[];

for i=1:size(StaInfo,2)
    n = n + 1;
    stalon=[stalon StaInfo(n).lon ];
    stalat=[stalat StaInfo(n).lat ];
end

lon0=src_grd(1);lat0=src_grd(2);
rdis=[];
for i=1:length(stalon)
    rdis=[rdis distance11(lat0,lon0,stalat(i),stalon(i),6371)];
end
[I,ind]=sort(rdis);


%% plot station recordings
%max_t=200;
figure('color',[1 1 1]);
  ha=tight_subplot(nline,nrow,[0.04 0.04],[0.04 0.04],[0.04 0.04]);
  time_file='time.dat';
  t_max=t(end)/60;
for kkk=1:min(nsta,nline*nrow)
    axes(ha(kkk));
     i=ind(kkk);
    %data
    plot(t(1:size(data_obs,2))/60,data_obs(i,:),'r','LineWidth',1.5);  hold on;
   
   % prediction
    plot(t(1:size(tsu_pred3,2)-2)/60,tsu_pred3(i,1:end-2),'k','LineWidth',1.5);  hold on;
    %text(t_max/2*0.5,-max(tsu_pred1(i,1:100))*0.8,StaInfo(i).name);
    text(t_max/2*0.5,-1,StaInfo(i).name);

    xlim([0 t(end)/60]);
    %ylim([min(min(data_obs(i,:)),-2) max(max(data_obs(i,:)),2) ] )
    
    % window for inversion
    bar_length=max(max(tsu_pred3(i,1:total_time))*1.5,0.5);
    plot([sta_start(i) sta_start(i)]/60, [-bar_length  bar_length], '--','color',[0.5 0.5 0.5])
    plot([sta_end(i) sta_end(i)]/60, [-bar_length  bar_length], '--','color',[0.5 0.5 0.5])
    if i==1
       legend('recorded','predicted','Start of window','End of window');
    end
end


%% Initial water elevation of a given iteration (iter_plot)
load japan_coastline.dat;
trench=load('trench_japan');
figure('color',[1 1 1]);
if n_slice>1
     src_total=src_1(1,:,:)-src_1(1,:,:); % initialization
     src_total=squeeze(src_total);
end

for i=1:n_slice % plot all slices

load([result_folder '/src_inversion_' num2str(iter_plot) '.mat']);

 nlon=size(src_1,3);nlat=size(src_1,2);
 if n_slice>1
 src_total=src_total+squeeze(src_1(i,:,:));
 end
       if n_slice>1
           subplot(1,n_slice+1,i);
       end
imagesc(linspace(srclonmn,srclonmx,nlon),linspace(srclatmn,srclatmx,nlat),squeeze(src_1(i,:,:)));
hold on;
%end

caxis([caxis_h(1),caxis_h(2)]);
set(gca,'ydir','normal')
xlim(xrange_plot);
ylim(yrange_plot);
box on;
hold on;
colorbar;

% addd contour line
high_end=ceil(max(max(squeeze(src_1(i,:,:)))));
low_end=floor(min(min(squeeze(src_1(i,:,:)))));
v=low_end:gap_contour:high_end;
[cs,h]=contour(linspace(srclonmn,srclonmx,nlon),linspace(srclatmn,srclatmx,nlat),squeeze(src_1(i,:,:)),v,'k','ShowText','on','LineWidth',1);
hold on;

z=clabel(cs, h, 'FontSize',14);
plot(src_grd(1), src_grd(2),'rp','MarkerSize',15,'MarkerFaceColor','r');hold on;
plot(trench(:,1),trench(:,2),'color',[1 1 1],'LineWidth',1.5);
colorbar;
hold on;
plot(stalon,stalat,'wv');
hold on;plot(japan_coastline(:,1),japan_coastline(:,2),'black');
title( ['Initial elevation at ' num2str(time_slice(i)) ' s of iteration ' num2str(iter_plot) '(m)'])
end

% plot the total source
if n_slice>1
subplot(1,n_slice+1,n_slice+1);
imagesc(linspace(srclonmn,srclonmx,nlon),linspace(srclatmn,srclatmx,nlat),src_total);
hold on;
caxis([caxis_h(1),caxis_h(2)]);
set(gca,'ydir','normal')
xlim(xrange_plot);
ylim(yrange_plot);
box on;
hold on;
colorbar;

high_end=ceil(max(max(src_total)));
low_end=floor(min(min(src_total)));

v=low_end:gap_contour:high_end;

[cs,h]=contour(linspace(srclonmn,srclonmx,nlon),linspace(srclatmn,srclatmx,nlat),src_total,v,'k','ShowText','on','LineWidth',1);

z=clabel(cs, h, 'FontSize',14);
plot(src_grd(1), src_grd(2),'rp','MarkerSize',15,'MarkerFaceColor','r');hold on;
plot(trench(:,1),trench(:,2),'color',[1 1 1],'LineWidth',1.5);
colorbar;
hold on;
plot(stalon,stalat,'wv');
hold on;plot(japan_coastline(:,1),japan_coastline(:,2),'black');
title( ['Stacked source of iteration ' num2str(iter_plot) '(m)'])

end


%% initial guessing
figure('color',[1 1 1]);
 for i=1:n_slice
    if n_slice>1
       subplot(1,n_slice,i);
    end
 [x, y, dat] = COMCOT_readBinaryDataSnapshot([result_folder '/z_1_' num2str((total_time-time_slice(i))/dt,'%05d') '_adjoint.dat']);
 %load ./inversion/debug_init_initial_result.mat
 imagesc(x,y,dat); hold on;
 plot(trench(:,1),trench(:,2),'color',[1 1 1],'LineWidth',1.5);
 plot(src_grd(1), src_grd(2),'rp','MarkerSize',15,'MarkerFaceColor','r');hold on;
 colorbar;
 hold on;
 plot(stalon,stalat,'kv');
 hold on;plot(japan_coastline(:,1),japan_coastline(:,2),'black');
 % caxis([caxis_h(1),caxis_h(2)]);
 title(['Initial model at ' num2str(time_slice(i)) ' s (m)'])
 set(gca,'ydir','normal')
 xlim(xrange_plot);
 ylim(yrange_plot);
 box on;
 end


 %%  plot misfit vs iteration
figure('color',[1 1 1]);
misfit_all=[];
 for i=1:nit
    load([result_folder '/data_pred_' num2str(i) '.mat']);
    res_it=0;
    for j=1:nsta
       res=tsu_pred(j,ceil(sta_start(j)/dt)+1:ceil(sta_end(j)/dt)-2)-data_obs(j,ceil(sta_start(j)/dt)+1:ceil(sta_end(j)/dt)-2);
       res_it=res_it+sqrt(mean(res.^2));
    end
    misfit_all=[misfit_all res_it];
 end
 misfit_all=misfit_all/nsta;
 plot(1: nit, misfit_all);
 grid on;
 box on;
 xlabel('Iteration');
 ylabel('misfit (m)')


  