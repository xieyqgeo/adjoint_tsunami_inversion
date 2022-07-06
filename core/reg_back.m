function src_back=reg_back(src,nsrcx1,nsrcy1)
[nsrcy,nsrcx]=size(src);
src_back=zeros(nsrcy1,nsrcx1);
nx=floor(nsrcx1/nsrcx);
ny=floor(nsrcy1/nsrcy);
for i=1:nsrcy
	for j=1:nsrcx
		src_temp = src(i,j) * ones(ny,nx);
		src_back(1+(i-1)*ny:i*ny,1+(j-1)*nx:j*nx) = src_temp;
	end
end
