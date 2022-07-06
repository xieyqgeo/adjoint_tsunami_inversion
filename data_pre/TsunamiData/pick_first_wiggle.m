load data_obs_ori.mat data_obs StaInfo dt t

% pick the first wiggle
% filter and pick up start time (and trim)
flag_use_all=0; % not completed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

data_obs2=data_obs;
t1=[]; t2=[];
for i=nsta:-1:1
    figure(1022);plot(t,data_obs2(i,:));xlabel('time (s)');
    set(gcf,'position',[252         570        1320         306]);
    [starttime,dummy] = ginput(1); 
    [firstpeak,dummy] = ginput(1);
    idx1 = round(starttime/dt); idx2=round(firstpeak/dt);
    fprintf('%s%s%s\n','start time is at the grid ',num2str(starttime),'s');
    fprintf('%s%s%s\n','start time is at the grid ',num2str(firstpeak),'s');
    prid = -(starttime-firstpeak);
    fprintf('%s%s%s\n','period is :',num2str(prid),'s');
    aaaa=input('if delete the record press d else press any key\n','s');
    if(strcmp(aaaa,'d'))
        data_obs(i,:)=[];
        StaInfo(i) = [];
    else
        t1=[t1;starttime];
        t2=[t2,firstpeak];
        % cut window between t1 and t2
        cutwin = zeros(1,length(t));
        tpr = tukeywin(50,0.5);
        tpr1 = tpr(1:14); 
        tpr2 = tpr(37:50);
        cutwin(idx1:idx2)=1; 
        cutwin(idx1-14:idx1-1)=0; %tpr1; 
        cutwin(idx2+1:idx2+14)=0; %tpr2;
        data_obs(i,:)=data_obs(i,:).*cutwin;
        
    end
end

figure(1003); hold on;
xlabel('time (min)'); title('trimmed tsunami recording')
for i=1:size(data_obs,1)
    plot(t/60,data_obs(i,:)/max(data_obs(i,:))+i,'b');
    text(0,i,StaInfo(i).name);
end

save data_obs_mod.mat data_obs StaInfo dt t t1 t2
    