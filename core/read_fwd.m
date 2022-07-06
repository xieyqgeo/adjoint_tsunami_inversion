function data=read_fwd(nt,dt,nsta) %,fl,fh,flag_filter)
ntt=ceil(nt/dt)+1;
data=zeros(nsta,ntt);
for i=1:nsta
	fn=sprintf('%s%04d%s','Station',i,'.dat');
	fid = fopen(fn, 'rb');
	StartTag = fread(fid, 1, 'integer*4');
	NDataLength = fread(fid, 1, 'integer*4');
	NFaults = fread(fid, 1, 'integer*4');
	ComputeGreen = fread(fid, 1, 'integer*4');
	EndTag = fread(fid, 1, 'integer*4');

	if(ComputeGreen == 1)
		    nDatCol = NFaults+1;
		else
			nDatCol = 2;
	end
	dat = zeros(NDataLength, nDatCol);
			
	StartTag = fread(fid, 1, 'integer*4');
	t = fread(fid, NDataLength, 'real*8');
	dat(:,1) = t(:);
	EndTag = fread(fid, 1, 'integer*4');
		
	for iCol = 1:nDatCol-1
	    StartTag = fread(fid, 1, 'integer*4');
	    h = fread(fid, NDataLength, 'real*8');
	    dat(:,1+iCol) = h(:);
	    EndTag = fread(fid, 1, 'integer*4');
	end
					
	fclose(fid);

	data(i,:)=dat(1:ntt,2);
end

% if flag_filter == 1
%   [BB,AA]=butter(4,fl/((1/dt)/2),'high');
%   [BB2,AA2]=butter(4,fh/((1/dt)/2),'low');
%   for i=1:nsta
%      data(i,:)=filtfilt(BB,AA,data(i,:));
%      data(i,:)=filtfilt(BB2,AA2,data(i,:));
%   end
% end
