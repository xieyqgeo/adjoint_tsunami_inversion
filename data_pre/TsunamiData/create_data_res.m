function create_data_res(data,dt)
nsta=size(data,1);
ndata=size(data,2);
fid=fopen('data_res_matlab.bin','wb');
fwrite(fid,nsta,'real*8');
fwrite(fid,ndata,'real*8');
fwrite(fid,dt,'real*8');
for ista=1:nsta
	h=data(ista,:);
	fwrite(fid, h, 'real*8');
end	
fclose(fid);
