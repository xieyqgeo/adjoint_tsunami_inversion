function  [src_1,ss,alpha,dk1,gk1]=run_time_reversal_resi_5slice(xr,yr,time_slice,pur,nsrcx,nsrcy,nsrcx1,nsrcy1,ccc,src1,flt2,gk,dk,data_res,time_step,src_grd,staweight,it,nproc,total_time,srclonmn,srclonmx,srclatmn,srclatmx,hmin,hmax)
% run time reversal, claculate residual, read and store snapshot, do the
% inversion???

% calculate backward modelling
cd ../data_res;
create_data_res(data_res,time_step);
aaaaa=system(['./convert_data_res']); 
if(aaaaa~=0) error('system wrong'); end
cd ../reversal;
create_comcot_input_reversal('../data_res/',src_grd,staweight);
create_comcot_input_pcomcot_reversal_5slice(1,total_time,time_step,time_slice);
display(['running pcomcot for iteration ',num2str(it),' reversal...']);
    
aaaaa=system(['mpirun -n ',num2str(nproc),' ./pcomcot > comcot_iter',num2str(it),'.log']); 
if(aaaaa~=0) error('system wrong'); end
% already runned for debug not run this
    
n_slice=length(time_slice);

for i=1:n_slice
    aaaaa=system(['cp _xcoordinate01.dat _ycoordinate01.dat z_01_0',num2str(total_time-time_slice(i),'%05d'),'.dat ../inversion']); 
    if(aaaaa~=0) error('system wrong'); end
end

adjoint_wfd=zeros(length(time_slice),nsrcy,nsrcx);
%adjoint_cut=zeros(length(time_slice),length(xr),length(yr));

   for i=1:n_slice 
	aaaaa=system(['cp z_01_0',num2str(total_time-time_slice(i),'%05d'),'.dat ../inversion/z_',num2str(it),'_',num2str(total_time-time_slice(i),'%05d'),'_adjoint.dat']); 
	if(aaaaa~=0) error('system wrong'); end
	cd ../inversion;
	snapshotnm=['z_',num2str(it),'_',num2str(total_time-time_slice(i),'%05d'),'_adjoint.dat'];
	[xx,yy,adjoint_wfd_i]=COMCOT_readBinaryDataSnapshot(snapshotnm);
        %adjoint_wfd(i,:,:)=adjoint_wfd_i;
	[xr,yr,adjoint_cut_i]=cutinit([srclonmn,srclonmx],[srclatmn,srclatmx],xx,yy,adjoint_wfd_i);
        adjoint_wfd(i,:,:)=adjoint_cut_i;
   end
   %???????????????????????????????????????????????/didn't change 
   aaaaa=system(['mv _InitialElevation_interp.dat neg_gradient_',num2str(it),'.dat']); 
   if(aaaaa~=0) error('system wrong'); end
	gk1=zeros(length(time_slice),nsrcy,nsrcx);
	for i=1:n_slice
    		gk1(i,:,:) = - reg_smoo(squeeze(adjoint_wfd(i,:,:)),nsrcx,nsrcy);
	end
	dk1 = conjugate_direction(it,gk1,gk,dk); %???????????????????????????????????????/didn't change
 	src_1 = src1;
    if ~isempty(flt2)
        for i=1:n_slice
    		src_1(i,:,:) = filter2(flt2,squeeze(src_1(i,:,:)));    % gaussian regularization
        end
    end   
    	cwin = 1 -  tukeywin(ccc*2,0.5);
    	cwin = cwin(ccc+1:end);

	%cut
	for k=1:n_slice
		for i=1:size(src_1,2)
  			tmp = squeeze(src_1(k,i,:));  
  			index0 = find(tmp==0);
  			if(~isempty(index0))   
    				index00 = index0(end);
    				tmp(1:index00) = 0;   
    				tmp(index00+1:index00+ccc) = tmp(index00+1:index00+ccc) .* cwin;
  			end
  			src_1(k,i,:) = tmp;
		end
	end
	%cut
	ss = zeros(length(time_slice),nsrcy,nsrcx);
	g_mean = (sum(dk1(:).*dk1(:)))^0.5;  %???????????????????????????????????//
	for i=1:n_slice
		ss_i = reg_smoo(squeeze(src_1(i,:,:)),nsrcx,nsrcy);
    		ss(i,:,:) = ss_i;
	end
        s_mean = (sum(ss(:).*ss(:)))^0.5;
	alpha=s_mean/g_mean*pur;
	display(['s_mean=',num2str(s_mean),' g_mean=',num2str(g_mean),' alpha=',num2str(alpha)]);
    	%back tracking to find the numerical step length %?????????????

    	if(it<10)
		f1=2;
	else
		f1=1;
    	end
    
	ss1 = ss + alpha * f1 * dk1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
	for k=1:n_slice
		src_1(i,:,:) = reg_back(squeeze(ss1(i,:,:)),nsrcx1,nsrcy1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
	end
	src_1(src_1<hmin)=hmin;src_1(src_1>hmax)=hmax; % restrict the water elevation range
