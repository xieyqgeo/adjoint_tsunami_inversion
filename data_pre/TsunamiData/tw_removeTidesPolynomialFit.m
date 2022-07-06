function tw_removeTidesPolynomialFit(filelistname, startt, endt, strPause)
% tw_removeTidesPolynomialFit('Stations.ctl')
% tw_removeTidesPolynomialFit('Stations.ctl', 1); pause

npol = 25; %order of polynomial fitting
removets = 0; removete = 20; %remove seismic signals which cause mismatch

f = fopen(filelistname,'r');
n = 0;
while ~feof(f)
    s = fgetl(f);
    n = n + 1;
    [s1,s] = strtok(s); StaInfo(n).lon = str2double(s1);
    [s1,s] = strtok(s); StaInfo(n).lat = str2double(s1);
    StaInfo(n).name = strtrim(s);
end
fclose(f);

startt0 = -2500; endt0 = 3500;
if(nargin <= 2)
    startt = startt0; endt = endt0;
end

IsPause = 0;
if(nargin == 4 || nargin == 2)
    IsPause = 1;
end

for i = 1:n
    filename = [StaInfo(i).name, '.txt'];
    fout = [StaInfo(i).name, '.dat'];
    data = load(filename);
    t = data(:,1); h = data(:,2);
    
    for j = 1:length(t)
        if(t(j)/60 >= startt)
            j1 = j; break
        end
    end
    for j = length(t):-1:1
        if(t(j)/60 <= endt)
            j2 = j; break
        end
    end
    t = t(j1:j2); h = h(j1:j2);
    p = polyfit(t, h, npol); %first fit
    
    for j = 1:length(t)
        if(t(j)/60 >= removets)
            j1 = j; break
        end
    end
    for j = length(t):-1:1
        if(t(j)/60 <= removete)
            j2 = j; break
        end
    end
    h2 = h; h2(j1:j2) = polyval(p,t(j1:j2)); %remove signal and do a second fit
    
    p = polyfit(t, h2, npol); hpoly = polyval(p,t); hres = h-hpoly;
    
    close all
    figure; hold on;
    plot(t/60, h, 'b'); %plot(t/60, h2, 'k');
    plot(t/60, hpoly, 'r');
    title(filename);
    xlabel('time (min)'); ylabel('water elevation (m)');
    set(gcf,'Units','normal','Position',[0.05,0.35,0.45,0.45]);
    
    figure; hold on
    plot(t/60, hres, 'b');
    xlim([0 max(t)]/60);
    title(filename);
    grid on
    xlabel('time (min)');
    ylabel('water elevation residual (m)');
    set(gcf,'Units','normal','Position',[0.55,0.35,0.45,0.45]);
    
    dlmwrite(fout,[t hres],'delimiter',' ');
    disp(fout);
    if(IsPause == 1 && i ~= n), input(''), end
end