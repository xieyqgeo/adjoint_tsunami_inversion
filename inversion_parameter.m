% ----------------------------- %
% set up inversion parameters   %
% for tsunami adjoint inversion %
% read by driver_comcot.m       %
% ----------------------------- %
clear all;

% The time snapshot for inversion
% If assuming all source emerge at the first second, set time_slice=0;
time_slice=[0 200]; 


% observation data
obs_path = "data_pre/data_obs.mat";
load(obs_path);

% number of cores
nproc = 2;

% number of iterations
nit   = 10; 

% epicenter [lon, lat]
src_grd=[142.861, 38.1035]; 

% simulation region
lonmn = 139; 
lonmx =  147; 
latmn = 33; 
latmx = 43;  

% tsunami source box range
srclonmn = 141; 
srclonmx = 145; 
srclatmn = 35; 
srclatmx = 40; 

% min and max water elevation
hmin = -20;
hmax =  20;

% filtering band for tsunami recordings
% filtering switch -0 no -1 yes
flag_filter = 1; % important for S-net to remove seismic signals
% low frequency
fl = 1/1000;
% high frequency
fh = 1/100;

source_scale = 1;   % scale factor for down-sampling the source region
     
% number of stations
nsta = size(data_obs,1); 

% station weighting
wt = ones(1,nsta);
%wt = [1, 1, 1, 1];

flt_len = 15;        % The length of filter length to smooth the result before start the next iteration
                     % If you don't want to smooth the result at the end of
                     % each iterations, set flt_len to 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all important parameters are set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% simulation time (in seconds)
total_time = t(end)+dt;   % will be overrided by the data total time 

% time step of calculation (in seconds)
time_step = dt;        % will be overrided by the data time interval

% ---------- parameter check ----------- %
% ---------- do not change   ----------- %
% checking parameters

if(time_step ~= dt_compute)
	fprintf('%s\n',"warning: data time step mismatch given time step, override with data time step");
	time_step = dt_compute;
end
if(total_time ~= max(t))
	fprintf('%s\n',"warning: data time span mismatch given total time, override with data time span");
        total_time = max(t)+dt;
end

% check if need interpolation
nlen = round(total_time/time_step);
if(size(data_obs,2)~=nlen)
	t_compute = 0:time_step:total_time;
	data_obs2 = data_obs;
	data_obs = zeros(nsta,nlen);
	for i=1:nsta
		data_obs(i,:) = interp1(t,data_obs2(i,:),t_compute);
	end
end

% save data_obs for inversion use
save data_obs_inv.mat data_obs time_step total_time

% ------- fixed parameters (do not change) -------- %

addpath('./core');
addpath('./core/yuqing');
fwd_dir  = './forward';
bkwd_dir = './reversal';
data_dir = './data_pre';
inv_dir  = './inversion';

d2g = 111194.9;

% how many time steps to save frame
time_save = total_time;

% station weighting 
staweight = wt;

% weighting for compute misfit
ws = wt;

