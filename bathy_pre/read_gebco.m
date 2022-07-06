clear all;
close all;

% This script is to cut a small piece of bathymetry for simulation, which cover all
% the stations.

% intended simulation region
lat01 = 43;   % north end
lat02 = 33;   % south end
lon01 = 139;   % west end
lon02 = 147.5;    % east end

% intended grid size ratio (default: 15 arc seconds)
gsr = 6; % if gsr=1, it means a grid size of 15 arc second. if gsr=2, it means 2*15 arc seoncd.

% read and cut GEBCO data
lat = ncread('GEBCO_2020.nc','lat');
lon = ncread('GEBCO_2020.nc','lon');
elev = ncread('GEBCO_2020.nc','elevation');
%I=find(lon<0);
%lon(I)=lon(I)+360;

% processing
if(lon01>lon02)
    idx1 = find(abs(lon-lon01)<0.003); 
    idx1 = idx1(end);
    idx2 = find(abs(lon-180)<0.003);
    idy1 = find(abs(lat-lat02)<0.003); %50
    idy1 = idy1(end);
    idy2 = find(abs(lat-lat01)<0.003); %30
    idy2 = idy2(1);
    latt = lat(idy1:idy2);
    lon1 = lon(idx1:idx2);
    bath1 = elev(idx1:idx2,idy1:idy2);

    idx1 = find(abs(lon--180)<0.003);
    idx2 = find(abs(lon-lon02)<0.003); %70
    idx2 = idx2(1);
    bath2 = elev(idx1:idx2,idy1:idy2);
    lon2 = lon(idx1:idx2);
    bath = [bath1;bath2];
    lonn = [lon1;lon2+360];
else
    idx1 = find(abs(lon-lon01)<0.003); 
    idx1 = idx1(end);
    idx2 = find(abs(lon-lon02)<0.003);
    idx2 = idx2(1);
    lonn = lon(idx1:idx2);
    idy1 = find(abs(lat-lat02)<0.003); %50
    idy1 = idy1(end);
    idy2 = find(abs(lat-lat01)<0.003); %30
    idy2 = idy2(1);
    latt = lat(idy1:idy2);
    bath = elev(idx1:idx2,idy1:idy2);
end

bath2 = bath(1:gsr:end,1:gsr:end);
lat2 = latt(1:gsr:end);
lon2 = lonn(1:gsr:end);
flt = fspecial('gaussian',9,4.5);
bath2 = filter2(flt,bath2);

figure
imagesc(lon2,lat2,transpose(bath2));
set(gca,'ydir','normal');


% % fill the small channels and remove small islands
% 
% for ii=1:size(bath3,1)
%     temp = bath3(ii,:);
%     for jj=11:length(temp)-10
%         if(temp(jj)>0)  % small islands
%             isa = temp(jj-10:jj+10)>0;
%         if(temp(jj)<0)
%             if(temp(jj-1)>0&&temp(jj
            
bath3=bath2;
% trim very small islands
tsize = 20;
istrim = ones(size(bath2,1)-tsize/2,size(bath2,2)-tsize/2);
for i=tsize/2+1:size(bath2,1)-tsize/2
    
    for j=1+tsize/2:size(bath2,2)-tsize/2
        if(i==2928&&j==2836)
            disp(2836)
        end
        temp = bath2(i-tsize/2:i+tsize/2,j-tsize/2:j+tsize/2);
        temp2 = sum(sum(temp>0));
        if(temp2/tsize>0.1)
            bath3(i,j) = sum(sum(temp))/tsize/tsize;
        end
    end
end

figure
imagesc(lon2,lat2,transpose(bath3));
set(gca,'ydir','normal');

fid = fopen('layer01.xyz','w');
for j=1:length(lat2)
    for i=1:length(lon2)
        fprintf(fid,'%f %f %f\n',lon2(i),lat2(j),-bath2(i,j));
    end
end

fclose(fid);

fid = fopen('_xcoordinate01.dat','w');

    for i=1:length(lon2)
        fprintf(fid,'%f ',lon2(i));
    end
fclose(fid);

fid = fopen('_ycoordinate01.dat','w');

    for i=1:length(lat2)
        fprintf(fid,'%f ',lat2(i));
    end
fclose(fid);