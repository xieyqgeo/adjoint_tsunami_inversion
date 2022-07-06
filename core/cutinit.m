function [xr,yr,TRIcut] = cutinit(xlims,ylims,x,y,snapshot)
%[X,Y] = meshgrid(x,y);
%xr = xlims(1):x(2)-x(1):xlims(2);
%yr = ylims(2):y(2)-y(1):ylims(1);
%[XR,YR] = meshgrid(xr,yr);
%TRIcut = interp2(X,Y,snapshot,XR,YR);
%dlmwrite('_xcoordinate01_interp.dat',xr');
%dlmwrite('_ycoordinate01_interp.dat',yr');
%dlmwrite('_InitialElevation_interp.dat',TRIcut,'delimiter',' ');
%dx=x(2)-x(1); dy=y(2)-y(1);
dx=(x(end)-x(1))/(length(x)-1); dy=(y(end)-y(1))/(length(y)-1);
grd_xlim1=round((xlims(1)-x(1))/dx); grd_xlim2=round((xlims(2)-x(1))/dx);   % W-E
grd_ylim1=round((ylims(1)-y(1))/dy); grd_ylim2=round((ylims(2)-y(1))/dy);   % N-S
xr=x(1)+dx*(grd_xlim1:grd_xlim2);
yr=y(1)+dy*(grd_ylim1:grd_ylim2);
TRIcut=snapshot(grd_ylim1+1:grd_ylim2+1,grd_xlim1+1:grd_xlim2+1);
dlmwrite('_xcoordinate01_interp.dat',xr');
dlmwrite('_ycoordinate01_interp.dat',yr');
dlmwrite('_InitialElevation_interp.dat',TRIcut,'delimiter',' ');
end
