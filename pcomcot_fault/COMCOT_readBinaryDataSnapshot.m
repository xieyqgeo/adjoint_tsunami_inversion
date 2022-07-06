function [x, y, dat] = COMCOT_readBinaryDataSnapshot(fn)
% [x, y, dat] = COMCOT_readBinaryDataSnapshot(fn)

lbslash = strfind(fn,'/');
if(isempty(lbslash))
    path = './';
    fn0 = fn;
else
    lbslashlast = max(lbslash);
    path = fn(1:lbslashlast);
    fn0 = fn(lbslashlast+1:end);
end

n = length(fn0);
ilayer = str2num(fn0(n-12:n-11));
if(fn0(1) == '_')
    ilayer = str2num(fn0(n-5:n-4));
end
if(isempty(ilayer) == 1)
    ilayer = 1;
end

x = load([path,'_xcoordinate',sprintf('%02d',ilayer),'.dat']); ncol = length(x);
y = load([path,'_ycoordinate',sprintf('%02d',ilayer),'.dat']); nrow = length(y);

if((fn0(1) == 'M') || (fn0(1) == 'm'))
    ncol = ncol - 1; % for M
elseif((fn0(1) == 'N') || (fn0(1) == 'n'))
    nrow = nrow - 1; % for N
end

% ncol = ncol - 1; % for M
% nrow = nrow - 1; % for N

fid = fopen(fn, 'rb');
StartTag = fread(fid, 1, 'integer*4');
ncol = fread(fid, 1, 'integer*4'); %NX = ncol
nrow = fread(fid, 1, 'integer*4'); %NY = nrow
EndTag = fread(fid, 1, 'integer*4');
dat = zeros(nrow, ncol);
for i = 1:nrow
    StartTag = fread(fid, 1, 'integer*4');
    dattmp = fread(fid, ncol, 'real*8');
    dat(i,:) = dattmp(:);
    EndTag = fread(fid, 1, 'integer*4');
end
fclose(fid);
