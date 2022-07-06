function [data_res,res0,tsu_pred]=run_forward(xx,yy,latmn,latmx,lonmn,lonmx,AA,BB,time_slice,data_obs_filter,aaa,it,src_1,nproc,total_timestep,time_step,nsta,srclonmn,srclonmx,srclatmn,srclatmx,ws,specweight)

total_time_i=[];
for i=1:length(time_slice)-1
  total_time_i=[total_time_i time_slice(i+1)-time_slice(i)];
end
total_time_i=[total_time_i total_timestep-time_slice(end)];


save_time_i=total_time_i;

last_wavefield=zeros(size(src_1,2),size(src_1,3));%%%%%
tsu_pred=[];
scale_y=(srclatmx-srclatmn)/size(src_1,2);
scale_x=(srclonmx-srclonmn)/size(src_1,3);
latindex1 = ceil((srclatmn-latmn)/scale_y);
latindex2 = ceil((srclatmx-latmn)/scale_y);
lonindex1=ceil((srclonmn-lonmn)/scale_x);
lonindex2=ceil((srclonmx-lonmn)/scale_x);
src_1_large=zeros(length(yy),length(xx));
init_TRI_i=zeros(length(yy),length(xx));
for i=1:length(time_slice)
    %create initial condition
   % create_init(squeeze(src_1(i,:,:))+last_wavefield); % add it into src_1
   src_1_temp = squeeze(src_1(i,:,:));   
   src_1_large(latindex1:latindex2-1,lonindex1:lonindex2-1) = src_1_temp;
    % create_init(src_1_large+init_TRI_i);
    % aaaaa=system(['./conv_init_waterelev']);
   filename_write='InitialElevation.xyz';
   fid=fopen(filename_write,'wt');%
    
   init_TRI_current=init_TRI_i+src_1_large;
   for jj=1:length(yy)
       for ii=1:length(xx)
            fprintf(fid,'%f %f %f\n',xx(ii),yy(jj),init_TRI_current(jj,ii)); 
       end 
   end
 
%    if(aaaaa~=0) error('system wrong'); end
    aaaaa=system(['cp InitialElevation.xyz ../forward']);	 
    if(aaaaa~=0) error('system wrong'); end
    cd ../forward;  
    %display(['save_time_i ',num2str(save_time_i)]);
    display(['i= ',num2str(i)]);
    create_comcot_input_pcomcot_forward_5slice(0,time_step,save_time_i(i));
	display(['running pcomcot for iteration ',num2str(it), ' ', num2str(i),'/' num2str(length(time_slice)),' forwarding...']);
	%aaaaa=system(['srun -n ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'_',num2str(i),'.log']); 
	aaaaa=system(['mpirun -n ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'_',num2str(i),'.log']);	
	if(aaaaa~=0) error('system wrong'); end
    aaaaa=system(['cp _xcoordinate01.dat _ycoordinate01.dat z_01_0',num2str(save_time_i(i),'%05d'),'.dat ../inversion']); 
    if(aaaaa~=0) error('system wrong'); end
    tsu_pred_i=read_fwd(total_time_i(i),time_step,nsta);   %%
    tsu_pred=[tsu_pred tsu_pred_i];
    cd ../inversion;
    snapshotnm=['z_01_0',num2str(save_time_i(i),'%05d'),'.dat'];
    [xx,yy,init_TRI_i]=COMCOT_readBinaryDataSnapshot(snapshotnm);
   % [xr,yr,last_wavefield]=cutinit([srclonmn,srclonmx],[srclatmn,srclatmx],xx,yy,init_TRI_i);
end

% read waveform
cd ../forward;

%calcualte residual
for kk=1:size(tsu_pred,1)
tsu_pred(kk,:)=filtfilt(BB,AA,squeeze(tsu_pred(kk,:)));
end	

data_res = (tsu_pred(:,1:size(data_obs_filter,2))-data_obs_filter).*repmat(ws',1,total_timestep).*specweight;
res0 = sum(data_res(:).*data_res(:));
display(['res0= ',num2str(res0)]);
