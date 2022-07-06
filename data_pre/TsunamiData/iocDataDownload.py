#! /usr/bin/python
#! original by Chao An, TZ revised for python3
import urllib, sys
import urllib.request

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


#stations = ['aric','pisa','iqui','pata','toco3','meji','mata','anto','papo','talt3',\
#'chnr2','sanf','juan','call','cald','huas2','coqu2','pich2','valp','sano','const',\
#'quir','talc','crnl','lebu'] #Iquique 2014 event
stations = ['kepo1','ohig','stan','imbt','salv','sthl2','mais','syow','kerg2'] 

#stations=['qtro','buca'] #no data!
StartTime = ['2021', '8', '12']
EndTime = ['2021', '8', '14']
EpiTime = ['2021','8','12','18','32','00']

#stations = ['lali','tala','balt','sant','tumc','buve2','malp2']  # Ecuador
#StartTime = ['2016', '4', '16']
#EndTime = ['2016', '4', '18']
#EpiTime = ['2016','4','16','23','58','36']
#stations = ['pich3','qtro2','valp','sano2','buca2','const','quir','talc2','crnl2','lebu2',\
#'coqu2','huas3','cald2','chnr2','talt3','papo2','anto2','meji2','toco3','pata2',\
#'juan2','sanf3'] #Southern Chile 2015
#StartTime = ['2015', '9', '16']
#EndTime = ['2015', '9', '18']
#EpiTime = ['2015','9','16','22','54','32']

tlen = -(secondsSince20000101(StartTime+['0','0','0'])-secondsSince20000101(EndTime+['0','0','0']))/24/3600.0
t0 = secondsSince20000101(EpiTime)

iocStationsTxt = ''

for IS in range(0,len(stations)):
	station = stations[IS]
	sys.stdout.write('Station '+station+':    ')

	urlioc = "http://www.ioc-sealevelmonitoring.org/bgraph.php?code=" + station + "&output=tab&period=" + str('%4.2f'%tlen) + "&endtime=" + '-'.join(EndTime)

	req = urllib.request.Request(urlioc)
	response = urllib.request.urlopen(req)
	the_page = response.read().lower().decode('utf-8')

	l = the_page.split('class=field>')
	i0 = -1
	for i in range(1,len(l)):
		if(l[i].split('<')[0].find('rad') != -1):
			i0 = i; break

	s = ''
	l = the_page.split('</td></tr><tr><td>')
	for i in range(1,len(l)):
		ll = l[i].split('</td><td>')
		h = ll[i0].strip()

		t = []
		ll = ll[0].strip().split('-')
		t.append(ll[0]); t.append(ll[1]);
		ll = ll[2].split(' ')
		t.append(ll[0])
		ll = ll[1].split(':')
		t = t+ll
		t = secondsSince20000101(t)
		try:
			float(h)
			s = s+str(t-t0)
			s = s+' '
			s = s+h+'\n'
		except ValueError:
			continue
	if(station[len(station)-1].isdigit()):
		stationName = station[0:len(station)-1]
	else:
		stationName = station
	f = open(stationName+'.txt','w')
	f.write(s)
	f.close()

	urlioc = "http://www.ioc-sealevelmonitoring.org/station.php?code="+station
	req = urllib.request.Request(urlioc)
	response = urllib.request.urlopen(req)
	the_page = response.read().lower().decode('utf-8')
	lon = the_page.split('longitude')[1].strip().split('</td></tr>')[0].split('class=nice>')[1]
	lat = the_page.split('latitude')[1].strip().split('</td></tr>')[0].split('class=nice>')[1]
	iocStationsTxt = iocStationsTxt + str('%5.3f'%(float(lon)))+'\t'+str('%5.3f'%(float(lat)))+'\t'+\
	stationName+'\n'

	sys.stdout.write('done\n')

f = open('iocStations.ctl','w')
f.write(iocStationsTxt)
f.close()
