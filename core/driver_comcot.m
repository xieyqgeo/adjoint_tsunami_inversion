% ------------------------------------------------------------------
% Least-squares RTM for tsunami initial water elevation inversion
% using COMCOT
% For information see Zhou et al., 2019. doi:10.1029/2018JB016678
% Contributors: Tong Zhou (tzhou@epss.ucla.edu) Yuqing Xie
% Last update: Sep 10, 2021
% ------------------------------------------------------------------
% set up parameters
inversion_parameter;
display(['Calculating para: dt=',num2str(time_step),';nt=',num2str(total_time)]);
display(['preparing data...']);

% window weight
specweight = ones(nsta,nlen);
for i=1:nsta
    specweight(i,1:floor(timewindow_start(i)/time_step)) = 0;
    specweight(i,ceil(timewindow_end(i)/time_step):end) = 0;
end

% first time reversal -> get initial guess
cd ./reversal;

create_comcot_input_reversal('../data_pre/',src_grd,staweight);
create_comcot_input_pcomcot(1,total_time,time_step,time_save);
display(['calculating time reversal as the initial guess...']);
display(['station weighting is\n',num2str(staweight)]);

display(['running pcomcot for 1st time reversal...']);
syscmd=system(['mpirun -np ',num2str(nproc),' ./pcomcot > comcot_reversal.log']); 
if(syscmd~=0) error('system wrong'); end

syscmd=system(['cp _xcoordinate01.dat _ycoordinate01.dat z_01_',num2str(ceil(total_time/time_step),'%06d'),'.dat ../inversion']); 
if(syscmd~=0) error('system wrong'); end

% set the initial guess
cd ../inversion;
snapshotnm=['z_01_',num2str(ceil(total_time/time_step),'%06d'),'.dat'];
[xx,yy,init_TRI]=COMCOT_readBinaryDataSnapshot(snapshotnm);
save debug_init_initial_result.mat xx yy init_TRI;
[xr,yr,init_guess]=cutinit([srclonmn,srclonmx],[srclatmn,srclatmx],xx,yy,init_TRI);
save degub_cut_initial_result.mat xr yr init_guess;
nsrcx1=length(xr); nsrcy1=length(yr);
nsrcx=length(xr); nsrcy=length(yr);% low resolution for inversion
save init_guess_2.mat init_guess;
%src_1=zeros(nsrcy,nsrcx);
%src_1(81:200,5:70)=1;
load('init_guess_2.mat'); %load initial model
ss1 = init_guess;
%ss1 = reg_smoo(init_guess,nsrcx,nsrcy);
src_1 = reg_back(ss1,nsrcx1,nsrcy1);
%src_1(src_1>0)=0;

% prepare for first forward propagation
cd ../forward;
create_comcot_input_pcomcot(0,total_time,time_step,time_save);
display(['prepare for iterations...\n']);
% iterations
misfit=zeros(nit+1,1);
src_nit=zeros(nsrcy1,nsrcx1,nit+1);
pur=0.05;
gk=zeros(nsrcy,nsrcx);
dk=zeros(nsrcy,nsrcx);
%Iteration Starts
for it=1:nit
	tic;
	if(it==1)
		thrsstep=1; pur=1;
	elseif(it<10&&it>1)
		thrsstep=0.1; pur=0.1;
	else
		thrsstep=0.05; pur=0.05;
	end
	% calculate forward modelling
	display(['Iteration ',num2str(it)]);
	cd ../inversion;
	create_init(src_1);
	syscmd=system(['./conv_init_waterelev']); 
	if(syscmd~=0) error('system wrong'); end
	syscmd=system(['cp InitialElevation.xyz ../forward']);	 
	if(syscmd~=0) error('system wrong'); end
	cd ../forward;
	display(['running pcomcot for iteration ',num2str(it),' forwarding...']);
	syscmd=system(['mpirun -np ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'.log']); 
	if(syscmd~=0) error('system wrong'); end
	tsu_pred=read_fwd(total_time,time_step,nsta,fl,fh,flag_filter);  
	data_res = (tsu_pred-data_obs).*repmat(ws', 1, nlen).*specweight;
	res0 = sum(data_res(:).*data_res(:));
	display(['res0= ',num2str(res0)]);
	% calculate backward modelling
	cd ../data_res;
	create_data_res(data_res,time_step);
	syscmd=system(['./convert_data_res']); 
	if(syscmd~=0) error('system wrong'); end
	cd ../reversal;
	create_comcot_input_reversal('../data_res/',src_grd,staweight);
	display(['running pcomcot for iteration ',num2str(it),' reversal...']);
	syscmd=system(['mpirun -np ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'.log']); 
	if(syscmd~=0) error('system wrong'); end
	syscmd=system(['cp z_01_',num2str(ceil(total_time/time_step),'%06d'),'.dat ../inversion/z_',num2str(it),'_adjoint.dat']); 
	if(syscmd~=0) error('system wrong'); end
	cd ../inversion;
	snapshotnm=['z_',num2str(it),'_adjoint.dat'];
	[xx,yy,adjoint_wfd]=COMCOT_readBinaryDataSnapshot(snapshotnm);
	[xr,yr,adjoint_cut]=cutinit([srclonmn,srclonmx],[srclatmn,srclatmx],xx,yy,adjoint_wfd);
	syscmd=system(['mv _InitialElevation_interp.dat neg_gradient_',num2str(it),'.dat']); 
	if(syscmd~=0) error('system wrong'); end
	save(['data_pred_',num2str(it),'.mat'],'tsu_pred');
	
	gk1 = - reg_smoo(adjoint_cut,nsrcx,nsrcy);
	dk1 = conjugate_direction(it,gk1,gk,dk);
	
	% set up the trial step length
	ss = reg_smoo(src_1,nsrcx,nsrcy);
	s_mean=(sum(ss(:).*ss(:)))^0.5;
	g_mean=(sum(dk1(:).*dk1(:)))^0.5;
	alpha=s_mean/g_mean*pur;
	display(['s_mean=',num2str(s_mean),' g_mean=',num2str(g_mean),' alpha=',num2str(alpha)]);
	
	%back tracking to find the numerical step length
	if(it<10)
		f1=2;
	else
		f1=1;
	end

	ss1 = ss + alpha * f1 * dk1;
	src_1 = reg_back(ss1,nsrcx1,nsrcy1);
	src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax; % restrict the water elevation range
	% second forward modelling
	create_init(src_1);
	syscmd=system(['./conv_init_waterelev']); 
	if(syscmd~=0) error('system wrong'); end
	syscmd=system(['cp InitialElevation.xyz ../forward']); 
	if(syscmd~=0) error('system wrong'); end
	cd ../forward;
	display(['running pcomcot for iteration ',num2str(it),' back tracking...']);
	syscmd=system(['mpirun -np ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'.log']); 
	if(syscmd~=0) error('system wrong'); end
	tsu_pred=read_fwd(total_time,time_step,nsta,fl,fh,flag_filter);   %%
	data_res = (tsu_pred - data_obs).*repmat(ws', 1, nlen).*specweight;
	res1 = sum(data_res(:).*data_res(:));
	display(['f1= ',num2str(f1),' res1= ',num2str(res1)]);
	if res1>res0
		while res1>res0 && f1>thrsstep
			f2=f1; res2=res1;
			f1 = f1*0.5;
			ss1 = ss + alpha * f1 * dk1;  % adjust searching step
			src_1 = reg_back(ss1,nsrcx1,nsrcy1);
			src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax;
			cd ../inversion;
			create_init(src_1);
			syscmd=system(['./conv_init_waterelev']); 
			if(syscmd~=0) error('system wrong'); end
			syscmd=system(['cp InitialElevation.xyz ../forward']); 
			if(syscmd~=0) error('system wrong'); end
			cd ../forward;
			display(['running pcomcot for iteration ',num2str(it),' back tracking...']);
			syscmd=system(['mpirun -np ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'.log']); 
			if(syscmd~=0) error('system wrong'); end
			tsu_pred=read_fwd(total_time,time_step,nsta,fl,fh,flag_filter);
			data_res = (tsu_pred - data_obs).*repmat(ws', 1, nlen).*specweight;
			res1 = sum(data_res(:).*data_res(:));
			display(['f1= ',num2str(f1),' res1= ',num2str(res1)]);
		end
	else
		f2=f1*2;
		ss1 = ss + alpha * f2 * dk1; % adjust searching step
		src_1 = reg_back(ss1,nsrcx1,nsrcy1);
		src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax;
		cd ../inversion;
		create_init(src_1);
		syscmd=system(['./conv_init_waterelev']); 
		if(syscmd~=0) error('system wrong'); end
		syscmd=system(['cp InitialElevation.xyz ../forward']); 
		if(syscmd~=0) error('system wrong'); end
		cd ../forward;
		display(['running pcomcot for iteration ',num2str(it),' back tracking...']);
		syscmd=system(['mpirun -np ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'.log']); 
		if(syscmd~=0) error('system wrong'); end
		tsu_pred=read_fwd(total_time,time_step,nsta,fl,fh,flag_filter);
		data_res = (tsu_pred - data_obs).*repmat(ws', 1, nlen).*specweight;
		res2 = sum(data_res(:).*data_res(:));
		display(['f2= ',num2str(f2),' res2= ',num2str(res2)]);
	end
	gama=(f1^2*(res0-res2)+f2^2*(res1-res0))/(2*res0*(f1-f2)+2*res1*f2-2*res2*f1);
	display(['gama= ',num2str(gama),' numerical step_length= ',num2str(gama*alpha)]);
	ss1 = ss + alpha * gama * dk1;   % actually update source with step gama and conjugate direction dk1
	src_1 = reg_back(ss1,nsrcx1,nsrcy1);
	src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax;
	cd ../inversion;
	create_init(src_1);
	syscmd=system(['./conv_init_waterelev']); 
	if(syscmd~=0) error('system wrong'); end
	syscmd=system(['cp InitialElevation.xyz ../forward']); 
	if(syscmd~=0) error('system wrong'); end
	cd ../forward;
	display(['running pcomcot for iteration ',num2str(it),' final forwarding...']);
	syscmd=system(['mpirun -np ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'.log']); 
	if(syscmd~=0) error('system wrong'); end
	tsu_pred=read_fwd(total_time,time_step,nsta,fl,fh,flag_filter);
	data_res = (tsu_pred - data_obs).*repmat(ws', 1, nlen).*specweight;
	res3 = sum(data_res(:).*data_res(:));
	display(['res3= ',num2str(res3)]);
	if (res3>res1 || res3>res2)  % updated source worse? take the better one of f1 and f2
		if res1>res2
			res0=res2;
			lamta=f2;
		else
			res0=res1;
			lamta=f1;
		end
		ss1 = ss + alpha * lamta * dk1;
		src_1 = reg_back(ss1,nsrcx1,nsrcy1);
		src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax;
	else
		res0=res3;
	end
	% Refresh
	gk = gk1;
	dk = dk1;
	cd ../inversion;
	save(['src_inversion_',num2str(it),'.mat'],'src_1');
	src_nit(:,:,it) = src_1;
	toc;
    misfit(nit+1)=res0;
end

cd ../inversion
save tsunami_inversion_result.mat
