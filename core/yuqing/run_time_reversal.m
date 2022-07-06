function [nsrcx,nsrcy,nsrcx1,nsrcy1]=run_time_reversal(time_save,src_grd, time_step,time_slice,staweight,total_time,nproc,srclonmn,srclonmx,srclatmn,srclatmx,source_scale)
% time reversal to get the first snapshot
create_comcot_input_reversal('../data_pre/',src_grd,staweight); % need to change %%%%%%%%%%%%%%
create_comcot_input_pcomcot(1,total_time,time_step,time_save); 
display(['calculating time reversal as the initial guess... \n']);
display(['station weighting is\n',num2str(staweight),'\n']);
display(['running pcomcot for 1st time reversal...']);
aaaaa=system(['mpirun -n ',num2str(nproc),' ./pcomcot > comcot_reversal.log']); 
if(aaaaa~=0) error('system wrong'); end

aaaaa=system(['cp _xcoordinate01.dat _ycoordinate01.dat z_01_0',num2str(total_time),'.dat ../inversion']); 
if(aaaaa~=0) error('system wrong'); end 

cd ../inversion;

snapshotnm=['z_01_0',num2str(total_time),'.dat'];
[xx,yy,init_TRI]=COMCOT_readBinaryDataSnapshot(snapshotnm);
save debug_init_initial_result.mat xx yy init_TRI;
[xr,yr,init_guess]=cutinit([srclonmn,srclonmx],[srclatmn,srclatmx],xx,yy,init_TRI);
save degub_cut_initial_result.mat xr yr init_guess;
nsrcx1=length(xr); nsrcy1=length(yr);
nsrcx=round(nsrcx1/source_scale); nsrcy=round(nsrcy1/source_scale); %% set desampling rate for source inversion
save init_guess_2.mat init_guess;
%%%%%%%%%%%%%%% need to change
