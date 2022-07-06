% ------------------------------------------------------------------
% Least-squares RTM for tsunami initial water elevation inversion
% using COMCOT
% For information see Zhou et al., 2019. doi:10.1029/2018JB016678
% Contributors: Tong Zhou (tzhou@epss.ucla.edu) Yuqing Xie
% Last update: June 26, 2022
% ------------------------------------------------------------------
% set up parameters

inversion_parameter;
display(['Calculating para: dt=',num2str(time_step),';nt=',num2str(total_time)]);
display(['preparing data...']);

%window weight
specweight = ones(nsta,nlen);
for i=1:nsta
    specweight(i,1:floor(timewindow_start(i)/time_step)) = 0;
    specweight(i,ceil(timewindow_end(i)/time_step):end) = 0;
end

%parameter initiation for time slice inversion
time_start_forward=0;
time_start_reversal=0;



%new parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
aaaaa=0; %%%%%%%%%%%%%%%%
if flt_len>0
   flt2=fspecial('gaussian',flt_len);  % generate gaussian filter template
else
   flt2=[];
end
total_timestep=round(total_time/time_step);


% first time reversal -> get initial guess
cd ./reversal;

create_comcot_input_reversal('../data_pre/',src_grd,staweight);

create_comcot_input_pcomcot_reversal_5slice(0,total_time,time_step,time_slice);

display(['calculating time reversal as the initial guess... \n']);
display(['station weighting is\n',num2str(staweight),'\n']);

display(['running pcomcot for 1st time reversal...']);
aaaaa=system(['mpirun -np ',num2str(nproc),' ./pcomcot > comcot_reversal.log']); 

if(aaaaa~=0) error('system wrong'); end
n_slice=length(time_slice);

for i=1:n_slice
aaaaa=system(['cp _xcoordinate01.dat _ycoordinate01.dat z_01_0',num2str(total_time-time_slice(i),'%05d'),'.dat ../inversion']); 
if(aaaaa~=0) error('system wrong'); end
end

% set the initial guess
cd ../inversion;

snapshotnm=['z_01_0',num2str(total_time-time_slice(1),'%05d'),'.dat'];
[xx,yy,init_TRI_i]=COMCOT_readBinaryDataSnapshot(snapshotnm);
[xr,yr,init_guess_i]=cutinit([srclonmn,srclonmx],[srclatmn,srclatmx],xx,yy,init_TRI_i);
init_TRI=zeros(n_slice,length(yy),length(xx));
init_guess=zeros(n_slice,length(yr),length(xr));
init_TRI(1,:,:)=init_TRI_i;
init_guess(1,:,:)=init_guess_i;

if n_slice>1
for i=2:n_slice
snapshotnm=['z_01_0',num2str(total_time-time_slice(i),'%05d'),'.dat'];
[xx,yy,init_TRI_i]=COMCOT_readBinaryDataSnapshot(snapshotnm);
init_TRI(i,:,:)=init_TRI_i;
[xr,yr,init_guess_i]=cutinit([srclonmn,srclonmx],[srclatmn,srclatmx],xx,yy,init_TRI_i);
init_guess(i,:,:)=init_guess_i;
%save degub_cut_initial_result.mat xr yr init_guess;
end
end

nsrcx1=length(xr); nsrcy1=length(yr);
nsrcx=round(nsrcx1/source_scale); nsrcy=round(nsrcy1/source_scale); %% set desampling rate for source inversion
save init_guess_3.mat init_guess;

%ss1 = init_guess;
if ~isempty(flt2)
   for i=1:n_slice
     init_guess(i,:,:) = filter2(flt2,squeeze(init_guess(i,:,:)));   % filtering on the initial guess
   end
end

ccc = 20;                                % attenuate the number of grids offshore. please change ccc if you changed the bathymetry
cwin = 1 -  tukeywin(ccc*2,0.5);
cwin = cwin(ccc+1:end);
for k=1:n_slice
for i=1:size(init_guess(k,:,:),2)
  tmp = squeeze(init_guess(k,i,:));
  index0 = find(tmp==0);
  if(~isempty(index0))
    index00 = index0(end);
    tmp(1:index00) = 0;
    tmp(index00+1:index00+ccc) = tmp(index00+1:index00+ccc) .* cwin;
  end
  init_guess(k,i,:) = tmp;
end
end
init_guess = init_guess * 2; % original is 10
src_1=init_guess;
for i=1:n_slice
ss1 = reg_smoo(squeeze(init_guess(i,:,:)),nsrcx,nsrcy);  %% down sampling the source region
src_1(i,:,:)=reg_back(ss1,nsrcx1,nsrcy1);
end
save current_result.mat; %(time reversal)
% load current_result.mat;


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

% initialization
misfit=zeros(nit+1,1);
src_nit=zeros(n_slice,nsrcy1,nsrcx1,nit+1);
pur=0.05;
gk=zeros(n_slice,nsrcy,nsrcx);
dk=zeros(n_slice,nsrcy,nsrcx);
fl=1/100;sr=1/time_step;
data_obs_filter=data_obs;
[BB,AA]=butter(4,fl/(sr/2),'low');
[BB2,AA2]=butter(4,fh/(sr/2),'high');
for kk=1:size(data_obs,1)
data_obs_filter(kk,:)=filtfilt(BB,AA,squeeze(data_obs(kk,:)));
data_obs_filter(kk,:)=filtfilt(BB2,AA2,data_obs_filter(kk,:));
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

%Iteration Starts
for it=1:nit
	tic;
	%if(it==1)
	%	aaa=1; pur=1;
	%elseif(it<10&&it>1)
	%	aaa=0.1; pur=0.1;
	%else
		aaa=0.05; pur=0.05;
        %end
    % calculate forward modelling
    display(['Iteration ',num2str(it)]);
	cd ../inversion;
    [data_res,res0,tsu_pred]=run_forward(xx,yy,latmn,latmx,lonmn,lonmx,AA,BB,time_slice,data_obs_filter,aaa,it,src_1,nproc,total_timestep,time_step,nsta,srclonmn,srclonmx,srclatmn,srclatmx,ws,specweight);
    cd ../inversion;
    save(['data_pred_',num2str(it),'.mat'],'tsu_pred'); 
    %load('../reversal/cur_state.mat');
    [src_1,ss,alpha,dk1,gk1]=run_time_reversal_resi_5slice(xr,yr,time_slice,pur,nsrcx,nsrcy,nsrcx1,nsrcy1,ccc,src_1,flt2,gk,dk,data_res,time_step,src_grd,staweight,it,nproc,total_time,srclonmn,srclonmx,srclatmn,srclatmx,hmin,hmax);
	
    % second forward modelling
    [data_res,res1,tsu_pred]=run_forward(xx,yy,latmn,latmx,lonmn,lonmx,AA,BB,time_slice,data_obs_filter,aaa,it,src_1,nproc,total_timestep,time_step,nsta,srclonmn,srclonmx,srclatmn,srclatmx,ws,specweight);
	if(it<10)
                f1=1;
        else
                f1=0.5;
        end
	if res1>res0
		while res1>res0 && f1>aaa 
			f2=f1; res2=res1;
			f1 = f1*0.5;
			ss1 = ss + alpha * f1 * dk1;  % adjust searching step
            for mi=1:n_slice
			   src_1(mi,:,:) = reg_back(squeeze(ss1(mi,:,:)),nsrcx1,nsrcy1);
            end
			src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax;
			cd ../inversion;
             [data_res,res1,tsu_pred]=run_forward(xx,yy,latmn,latmx,lonmn,lonmx,AA,BB,time_slice,data_obs_filter,aaa,it,src_1,nproc,total_timestep,time_step,nsta,srclonmn,srclonmx,srclatmn,srclatmx,ws,specweight);
		end
    else
        
		f2=f1*2;
		ss1 = ss + alpha * f2 * dk1; % adjust searching step
        for mi=1:n_slice
		  src_1(mi,:,:) = reg_back(squeeze(ss1(mi,:,:)),nsrcx1,nsrcy1);
        end
		src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax;
		cd ../inversion;
         [data_res,res2,tsu_pred]=run_forward(xx,yy,latmn,latmx,lonmn,lonmx,AA,BB,time_slice,data_obs_filter,aaa,it,src_1,nproc,total_timestep,time_step,nsta,srclonmn,srclonmx,srclatmn,srclatmx,ws,specweight);
    end
    
	gama=(f1^2*(res0-res2)+f2^2*(res1-res0))/(2*res0*(f1-f2)+2*res1*f2-2*res2*f1);
	display(['gama= ',num2str(gama),' numerical step_length= ',num2str(gama*alpha)]);
	ss1 = ss + alpha * gama * dk1;   % actually update source with step gama and conjugate direction dk1
    for mi=1:n_slice
	   src_1(mi,:,:) = reg_back(squeeze(ss1(mi,:,:)),nsrcx1,nsrcy1);
    end

	src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax;
	cd ../inversion;
 [data_res,res3,tsu_pred]=run_forward(xx,yy,latmn,latmx,lonmn,lonmx,AA,BB,time_slice,data_obs_filter,aaa,it,src_1,nproc,total_timestep,time_step,nsta,srclonmn,srclonmx,srclatmn,srclatmx,ws,specweight);
	if (res3>res1 || res3>res2)  % updated source worse? take the better one of f1 and f2
		if res1>res2
			res0=res2;
			lamta=f2;
		else
			res0=res1;
			lamta=f1;
		end
		ss1 = ss + alpha * lamta * dk1;
    for mi=1:n_slice
	   src_1(mi,:,:) = reg_back(squeeze(ss1(mi,:,:)),nsrcx1,nsrcy1);
    end
    src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax;
	else
		res0=res3;
	end
	% Refresh
	gk = gk1;
	dk = dk1;
	cd ../inversion;
	save(['src_inversion_',num2str(it),'.mat'],'src_1');
    
	src_nit(:,:,:,it) = src_1;
	toc;
    misfit(nit+1)=res0;
    save(['record_iteration',num2str(it),'.mat']);
    end

save tsunami_inversion_result.mat
