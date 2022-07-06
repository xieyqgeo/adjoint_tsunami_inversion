#! /usr/bin/python

import urllib, urllib.request, sys, re, math, os

def secondsSince20000101(time0):
        for i in range(0,len(time0)):
                time0[i]=int(time0[i])
        Year=time0[0]; Month=time0[1]; Day=time0[2];
        Days = 0
        for i in range(2000, Year):
                Days += 365
                if(((i%4==0)&(i%100!=0)) | (i%400==0)):
                        Days += 1
        DaysInMonth=[31,28,31,30,31,30,31,31,30,31,30,31]
        if(((Year%4==0)&(Year%100!=0)) | (Year%400==0)):
                DaysInMonth[1] += 1
        for i in range(1,Month):
                Days += DaysInMonth[i-1]
        Days += Day-1
        return Days*24*3600+time0[3]*3600+time0[4]*60+time0[5]


stations = ['32401','32402','32412','32413']
StartTime = ['2014', '4', '1']
EndTime = ['2014', '4', '3']
time0 = [2014,4,1,23,46,47]





SpacePattern = re.compile('\s+')
if(len(stations) == 0):
        f = open('DARTStationsInfo.txt', 'r')
        while True:
                s = f.readline()
                if(len(s) == 0):
                        break
                stations.append(SpacePattern.split(s.strip())[0])
        f.close()

StationsGood = []
StationsTxt = ''
StartStr = "#yy  mm dd hh mm ss t   height"
EndStr = "</textarea></label></pre>"
t0 = secondsSince20000101(time0)

for IS in range(0,len(stations)):
        DataFileName = stations[IS]+'.txt'
        try:
                os.remove(DataFileName)
        except OSError:
                pass
        DataFileName = stations[IS]+'.dat'
        try:
                os.remove(DataFileName)
        except OSError:
                pass

for IS in range(0,len(stations)):
        station = stations[IS]
        sys.stdout.write('Station '+station+':    ')

        urlNOAADart = "http://www.ndbc.noaa.gov/station_page.php?station=" + station + "&type=0&startyear=" + StartTime[0] + "&startmonth=" + StartTime[1] + "&startday=" + StartTime[2] + "&endyear=" + EndTime[0] + "&endmonth=" + EndTime[1] + "&endday=" + EndTime[2] + "&submit=Submit"

        req = urllib.request.Request(urlNOAADart)
        response = urllib.request.urlopen(req)
        the_page = response.read().lower()

        if(the_page.find(StartStr) == -1):
                sys.stdout.write('NO DATA\n')
        else:
                s = re.compile('<b>\d+\.\d+\s+[ns]\s+\d+\.\d+\s+[ew]').findall(the_page)[0].replace('<b>','')
                if(s.find('n') != -1):
                        stnlat = s.split('n')[0].strip()
                else:
                        stnlat = '-'+s.split('s')[0].strip()
                if(s.find('e') != -1):
                        stnlon = s.replace('e','').strip().split(' ')[-1]
                else:
                        stnlon = '-'+s.replace('w','').strip().split(' ')[-1]

                s = StartStr +  the_page.split(StartStr)[1].split(EndStr)[0]
                l = s.split('\n')
                IsEmpty = True
                for i in range(2,len(l)):
                        ss = l[i].strip()
                        if((len(ss) > 0) and (ss.endswith('9999.000') == False)):
                                IsEmpty = False; break
                if(IsEmpty == True):
                        sys.stdout.write('NO DATA\n')
                else:
                        SourceFileName = stations[IS]+'.txt'
                        f = open(SourceFileName,'w')
                        f.write(s)
                        f.close()

                        SourceFile = open(SourceFileName,'r')
                        time = []; data = []
                        while True:
                                s = SourceFile.readline()
                                if len(s) == 0:
                                        break
                                s = s.strip()
                                if len(s) == 0:
                                        continue
                                if s.startswith('#'):
                                        continue
                                l = SpacePattern.split(s)
                                time1 = l[0:6]
                                t1 = secondsSince20000101(time1)
                                data1 = float(l[7])
                                if(math.fabs(data1) > 9000.0):
                                        continue
                                if(len(time)>0) and (t1-t0) == time[len(time)-1]:
                                        continue
                                if(len(data)>0 and math.fabs(data1-data[len(data)-1]) > 10.0):
                                        continue
                                time.insert(len(time),t1-t0)
                                data.insert(len(data),data1)
                        SourceFile.close()
                        DataFileName = stations[IS]+'.txt'
                        DataFile = open(DataFileName,'w')
                        for i in range(len(time)-1,-1,-1):
                                DataFile.write(str(time[i]).ljust(15,' ')+str(data[i])+'\n')
                        DataFile.close()

                        sys.stdout.write('done\n')
                        StationsGood.append(stations[IS])
                        StationsTxt = StationsTxt + stnlon+' '+stnlat+' '+stations[IS]+'\n'

f = open('DARTStations.ctl','w')
f.write(StationsTxt)
f.close()

