function src_smoo=reg_smoo(src,nsrcx,nsrcy)
[nsrcy1,nsrcx1]=size(src);
src_smoo=zeros(nsrcy,nsrcx);
nx=floor(nsrcx1/nsrcx);
ny=floor(nsrcy1/nsrcy);
for i=1:nsrcy
	for j=1:nsrcx
		src_temp = src(1+(i-1)*ny:i*ny,1+(j-1)*nx:j*nx);
		src_smoo(i,j) = sum(src_temp(:))/length(src_temp(:));
	end
end
