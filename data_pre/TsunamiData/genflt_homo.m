function outflt=genflt_homo(dt,fltlen,fcut)
%
    fnyq=1/dt/2;
    fper=fcut./fnyq;
    f=[0 fper fper+0.07 fper+0.14 1];
    m=[1 1 0.5 0 0];
    wd=kaiser(fltlen,6.5);
    bb=fir2(fltlen-1,f,m,wd);
    %bb2=ftrans2(bb);
    outflt=bb;
    
end