function init_guess=forward_simulation_first_prepare(flt2, ampli,ccc,nsrcx,nsrcy,total_time,time_step,time_save,nit,nsrcy1,nsrcx1,fl,sr)
load('init_guess_2.mat'); %load initial model

init_guess = filter2(flt2,init_guess);  % filtering on the initial guess

% attenuate the number of grids offshore. please change ccc if you changed the bathymetry
cwin = 1 -  tukeywin(ccc*2,0.5);
cwin = cwin(ccc+1:end);
for i=1:size(init_guess,1)
  tmp = init_guess(i,:);
  index0 = find(tmp==0);
  if(~isempty(index0))
    index00 = index0(end);
    tmp(1:index00) = 0;
    disp(index00);
    disp(ccc);
    disp(size(tmp));
    disp(size(init_guess))
    tmp(index00+1:index00+ccc) = tmp(index00+1:index00+ccc) .* cwin';
  end
  init_guess(i,:) = tmp;
end
init_guess = init_guess * ampli;
