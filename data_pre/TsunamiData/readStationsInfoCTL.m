function [StaInfo, nStations] = readStationsInfoCTL(fn)
% [StaInfo, nStations] = readStationsInfoCTL(fn)
% input: fn
% output: StaInfo(nStations), nStations
% StaInfo(i).name, StaInfo(i).lon, StaInfo(i).lat

f = fopen(fn, 'r');
n = 0;
while ~feof(f)
    s = fgetl(f);
    n = n + 1;
    [s1,s] = strtok(s); StaInfo(n).lon = str2double(s1);
    [s1,s] = strtok(s); StaInfo(n).lat = str2double(s1);
    StaInfo(n).name = strtrim(s);
end
fclose(f);
nStations = n;
