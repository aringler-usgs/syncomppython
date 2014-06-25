#!/usr/bin/env python

import sys
import os
import glob
import numpy
import matplotlib
import math
import warnings
from matplotlib.pyplot import (figure,axes,plot,xlabel,ylabel,title,subplot,legend,savefig,show)
from obspy import read, Stream
from obspy.core import UTCDateTime
from obspy.xseed import Parser
from time import gmtime, strftime
from obspy.core.util.geodetics import gps2DistAzimuth
from obspy.signal.cross_correlation import xcorr
from obspy.core.util.version  import read_release_version


debug=False
datalessloc = '/APPS/metadata/SEED/'
#Here is the data location use True for xs0 otherwise use false
dataloc = False
userminfre = .0025
usermaxfre = .01
lents = 4000
#Use half the value you think you want e.g. 2 gives you a total of 4 poles
filtercornerpoles = 2

newVer = True
#ver = (read_release_version()).split('-')[0]
#ver = ver.split('.')
#if int(ver[1]) > 8:
#	print 'Using new taper flag'
#	newVer = True

manstalist=False
stations=['IC BJT']

def getorientation(net,sta,loc,chan,evetime,xseedval):
#A function to get the orientation of a station at a specific time
	for cursta in xseedval.stations:
#As we scan through blockettes we need to find blockettes 50 and 52
		for blkt in cursta:
			if blkt.id == 50:
#Pull the station info for blockette 50
				stacall = blkt.station_call_letters.strip()
			if stacall == sta:
				if blkt.id == 52 and blkt.location_identifier == loc and blkt.channel_identifier == chan:
					if type(blkt.end_date) is str:
						curdoy = strftime("%j",gmtime())
						curyear = strftime("%Y",gmtime())
						curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
						if blkt.start_date <= evetime:
							azimuth = blkt.azimuth
					elif blkt.start_date <= evetime and blkt.end_date >= evetime:
						azimuth = blkt.azimuth
	return azimuth

def getdip(net,sta,loc,chan,evetime,xseedval):
#A function to get the dip of a station at a specific time
	for cursta in xseedval.stations:
#As we scan through blockettes we need to find blockettes 50 and 52
		for blkt in cursta:
			if blkt.id == 50:
#Pull the station info for blockette 50
				stacall = blkt.station_call_letters.strip()
			if stacall == sta:
				if blkt.id == 52 and blkt.location_identifier == loc and blkt.channel_identifier == chan:
					if type(blkt.end_date) is str:
						curdoy = strftime("%j",gmtime())
						curyear = strftime("%Y",gmtime())
						curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
						if blkt.start_date <= evetime:
							dip = blkt.dip
					elif blkt.start_date <= evetime and blkt.end_date >= evetime:
						dip = blkt.dip
	return dip

def rotatehorizontal(stream, angle):
	if stream[0].stats.channel in set(['LHE','LHN']):
		stream.sort(['channel'],reverse=False)
		stream[1].data = -stream[1].data
#Function to rotate the data by a given angle
        theta_r = math.radians(angle)
# create new trace objects with same info as previous
        rotatedN = stream[0].copy()
        rotatedE = stream[1].copy()
# assign rotated data
        rotatedN.data = stream[0].data*math.cos(- theta_r) + stream[1].data*math.sin(- theta_r)
        rotatedE.data = stream[1].data*math.cos(- theta_r) - stream[0].data*math.sin(- theta_r)
	rotatedN.stats.channel='LHN'
	rotatedE.stats.channel='LHE'
# return new streams object with rotated traces
        streamsR = Stream(traces = [rotatedN, rotatedE])
        return streamsR

def choptocommon(stream):
#A function to chop the data to a common time window
	stimes = []
	etimes = []

	for trace in stream:
		stimes.append(trace.stats.starttime)
		etimes.append(trace.stats.endtime)
	newstime = stimes[0]
	newetime = etimes[0]

	for curstime in stimes:
		if debug:
			print(curstime)
		if curstime >= newstime:
			newstime = curstime
	
	for curetime in etimes:
		if debug:		
			print(curetime)
		if curetime <= newetime:
			newetime = curetime

	if debug:
		print(newstime)
		print(newetime)
		print(stream)
	for trace in stream:	
		trace.trim(starttime=newstime,endtime=newetime)
	if debug:
		print(stream)
	return stream

def getlatlon(sta,etime,xseedval):
#A function to get the lat and lon of a station at a given time
	for cursta in xseedval.stations:
#As we scan through blockettes we need to find blockettes 50 and 52
		for blkt in cursta:
			if blkt.id == 50:
#Pull the station info for blockette 50
				stacall = blkt.station_call_letters.strip()
				if stacall == sta:
					lat = blkt.latitude
					lon = blkt.longitude	
					if type(blkt.end_effective_date) is str:
						curdoy = strftime("%j",gmtime())
						curyear = strftime("%Y",gmtime())
						curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
						if blkt.start_effective_date <= etime:
							lat = blkt.latitude
							lon = blkt.longitude
					elif blkt.start_effective_date <= etime and blkt.end_effective_date >= etime:
						lat = blkt.latitude
						lon = blkt.longitude	
	return lat,lon

def getstalist(sp,etime,curnet):
#A function to get a station list
	stations = []
	for cursta in sp.stations:
#As we scan through blockettes we need to find blockettes 50 
		for blkt in cursta:
			if blkt.id == 50:
#Pull the station info for blockette 50
				stacall = blkt.station_call_letters.strip()
				if debug:
					print "Here is a station in the dataless" + stacall
				if type(blkt.end_effective_date) is str:
					curdoy = strftime("%j",gmtime())
					curyear = strftime("%Y",gmtime())
					curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
					
					if blkt.start_effective_date <= etime:
						stations.append(curnet + ' ' + blkt.station_call_letters.strip())
				elif blkt.start_effective_date <= etime and blkt.end_effective_date >= etime:
					stations.append(curnet + ' ' + \
					blkt.station_call_letters.strip())	
	return stations

def readcmt(cmt):
	debugreadcmt = False
#Now we can continue like there is no difference between Princeton and our Synthetics
#Lets get the event time from the cmt
	cmtline1 = ' '.join(cmt[0].split())
	cmtlat = cmt[4].replace('latitude:','').strip()
	cmtlon = cmt[5].replace('longitude:','').strip()
	tshift = float(cmt[2].replace('time shift:','').strip())
	hdur = float(cmt[3].replace('half duration:','').strip())
	if debugreadcmt:
		print cmtline1
	cmtline1 = cmtline1.split()
	if debugreadcmt:
		print cmtline1[1] + ' ' + cmtline1[2] + ' ' + cmtline1[3] + ' ' + cmtline1[4] + ' ' + cmtline1[5] + ' ' + cmtline1[6]
	eventtime = UTCDateTime(int(cmtline1[1]),int(cmtline1[2]),int(cmtline1[3]),int(cmtline1[4]),int(cmtline1[5]),float(cmtline1[6]))
	if debugreadcmt:
		print 'Year:' + str(eventtime.year)
		print 'Day:' + str(eventtime.julday)
		print 'Hour:' + str(eventtime.hour)
		print 'Minute:' + str(eventtime.minute)
	return cmtlat, cmtlon, eventtime, tshift,hdur

def getdata(net,sta,eventtime,lents,dataloc):
#This function goes to one of the archives and gets the data
	debuggetdata = True
	preeventday = eventtime - 24*60*60
	posteventday = eventtime + 24*60*60
	prepostwin= 3000
#If II get off of /tr1 else get the data from /xs0 or /xs1
	if net == 'II':
		dataprefix1 = '/tr1/telemetry_days/'
	else:
		if net in set(['IW','NE','US']):	
			dataprefix = 'xs1'	
		else:
			dataprefix = 'xs0'
		dataprefix1 = '/' + dataprefix + '/seed/'
#	if not dataloc:
		dataprefix2 = '/tr1/telemetry_days/'
	if debug:
		print 'Here is the dataprefix:' + dataprefix
	try:
		st = read(dataprefix1 + net + '_' + sta + '/' + str(eventtime.year) + \
		'/' + str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '*/*LH*.seed', \
		starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
		st += read(dataprefix1 + net + '_' + sta + '/' + str(posteventday.year) + \
		'/' + str(posteventday.year) + '_' + str(posteventday.julday).zfill(3) + '*/*LH*.seed', \
		starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
		st += read(dataprefix1 + net + '_' + sta + '/' + str(preeventday.year) + \
		'/' + str(preeventday.year) + '_' + str(preeventday.julday).zfill(3) + '*/*LH*.seed', \
		starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
		
	except:
		st = read(dataprefix2 + net + '_' + sta + '/' + str(eventtime.year) + \
		'/' + str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '*/*LH*.seed', \
		starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
		st += read(dataprefix2 + net + '_' + sta + '/' + str(posteventday.year) + \
		'/' + str(posteventday.year) + '_' + str(posteventday.julday).zfill(3) + '*/*LH*.seed', \
		starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
		st += read(dataprefix2 + net + '_' + sta + '/' + str(preeventday.year) + \
		'/' + str(preeventday.year) + '_' + str(preeventday.julday).zfill(3) + '*/*LH*.seed', \
		starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
	st.merge(fill_value='latest')
	if debuggetdata:
		print 'We have data'
	return st

def getPAZ2(sp,net,sta,loc,chan,eventtime):
	debuggetPAZ2 = False
        data = {}
	station_flag = False
	channel_flag = False
	for statemp in sp.stations:
		for blockette in statemp:
			if blockette.id == 50:
				station_flag = False
				if net == blockette.network_code and sta == blockette.station_call_letters:
					station_flag = True
					if debuggetPAZ2:
						print 'We found the station blockettes'
			elif blockette.id == 52 and station_flag:
				channel_flag = False
				if blockette.location_identifier == loc and blockette.channel_identifier == chan:
					if debuggetPAZ2:
						print 'We are in the location and channel blockette'
						print 'End date: ' + str(blockette.end_date)
						print 'Start date: ' + str(blockette.start_date)
					if type(blockette.end_date) is str:
						curdoy = strftime("%j",gmtime())
						curyear = strftime("%Y",gmtime())
						curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
						if blockette.start_date <= eventtime:
							channel_flag = True
							if debuggetPAZ2:
								print 'We found the channel blockette'
					elif blockette.start_date <= eventtime and blockette.end_date >= eventtime:
						channel_flag = True
						if debuggetPAZ2:
							print 'We found the channel blockette'
			elif blockette.id == 58 and channel_flag and station_flag:
				if blockette.stage_sequence_number == 0:
					data['sensitivity'] = blockette.sensitivity_gain
				elif blockette.stage_sequence_number == 1:
					data['seismometer_gain'] = blockette.sensitivity_gain
				elif blockette.stage_sequence_number == 2:
					data['digitizer_gain'] = blockette.sensitivity_gain
			elif blockette.id == 53 and channel_flag and station_flag:
				data['gain'] = blockette.A0_normalization_factor
				data['poles'] = []
				if not blockette.transfer_function_types == 'A':
					msg = 'Only supporting Laplace transform response ' + \
					'type. Skipping other response information.'
					warnings.warn(msg, UserWarning)
					continue
				for i in range(blockette.number_of_complex_poles):
					p = complex(blockette.real_pole[i], blockette.imaginary_pole[i])
					data['poles'].append(p)
				data['zeros'] = []
				for i in range(blockette.number_of_complex_zeros):
					z = complex(blockette.real_zero[i], blockette.imaginary_zero[i])
					data['zeros'].append(z)
        return data

def getcolor(chan,loc):
	if chan in set(['LXN','LXE','LXZ']):
		color = 'k'
	elif (loc == '00' or loc ==''):
		color = 'g'
	elif loc == '10':
		color = 'r'
	elif loc == '60':
		color = 'c'
	else:
		color = 'b'
	return color

#def stfconv(syntrace,hdur):
	#A function to convolve the source time function with
#	stf = numpy.zeros(syntrace.stats.npts)
#	for ind in range(1,syntrace.stats.npts):
#		if (ind - 2*hdur) <= 0:
#			stf(ind) = 1
#	stf = stf*(1/(2*hdur))
	





#	return syntrace










def writestats(statfile,streamin,comp):
	try:
		syncomp = "LX" + comp	
		datacomp = "LH" + comp
		syn = streamin.select(channel = syncomp)
		for tr in streamin.select(channel = datacomp):	
			resi = "{0:.2f}".format(numpy.sum(tr.data*syn[0].data)/numpy.sum(numpy.square(syn[0].data)))
			lag, corr = xcorr(tr,syn[0],500)
			corr = "{0:.2f}".format(corr)
			statfile.write(tr.stats.network + "," + tr.stats.station)
			statfile.write("," + tr.stats.location + "," + tr.stats.channel + "," +  str(resi))
			statfile.write("," + str(lag) + "," + str(corr) + "\n")
	
	except:	
		if debug:
			print 'No residual for' + cursta + ' ' + 'LH' + comp	
	return




#Start of the main part of the program
if not len(sys.argv) == 4:
	print "Usage: Syntheticlocation ResultsName Network"
	exit(0)
 
synfile = sys.argv[1]
if synfile[-1] == '/':
	synfile = synfile[:-1]

#Read in the CMT solution 
if debug:
	print "We are using local synthetics"
if not os.path.isfile(synfile + '/CMTSOLUTION'):
	print "No CMT found"
	exit(0)
cmt = tuple(open(synfile + '/CMTSOLUTION'))
cmtlat, cmtlon, eventtime, tshift, hdur = readcmt(cmt)

#Lets make a local results directory
resultdir = sys.argv[2]
if resultdir[-1] == '/':
	resultdir = resultdir[:-1]

if not os.path.exists(os.getcwd() + '/' + resultdir):
	os.mkdir(os.getcwd() + '/' + resultdir)
evename = synfile.split("/")
evename = evename[len(evename)-1]

curnet = sys.argv[3]

statfile = open(os.getcwd() + '/' + resultdir + '/Results' + evename + curnet + '.csv' ,'w')
statfile.write('net,sta,loc,chan,scalefac,lag,corr\n')




#Lets read in the dataless
try:
	sp = Parser(datalessloc + curnet + ".dataless")
except:
	print "Can not read the dataless."
	exit(0)
if not manstalist:
	stations = getstalist(sp,eventtime,curnet)

if debug:
	print "Here are the stations we found"	
	for sta in stations:
		print "Here is a station:" + sta


#Lets start by using a station list and then move to a different approach
for sta in stations:
	cursta = sta.strip()
	if debug:
		print 'Current station:' + cursta
	cursta = sta.split()
	net = cursta[0]
	cursta = cursta[1]
	try:
		st = getdata(net,cursta,eventtime,lents,dataloc)
	except:
		print('No data for ' + net + ' ' + cursta)
		continue
		
#Lets go through each trace in the stream and deconvolve and filter
	for trace in st:
#Here we get the response and remove it
#Here is where I am mucking around		
		paz=getPAZ2(sp,net,cursta,trace.stats.location,trace.stats.channel,eventtime)
		try:
			if newVer:
				trace.taper(max_percentage=0.05, type='cosine')
			else:
				trace.taper()
			trace.simulate(paz_remove=paz)
#Here we filter
			trace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
			trace.integrate()
			if newVer:
				trace.taper(max_percentage=0.05, type='cosine')
			else:
				trace.taper()
			trace.trim(starttime=eventtime + tshift/2,endtime=(eventtime+lents + tshift/2))
			trace.detrend()
			trace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
		except:
			print('Can not find the response')
			st.remove(trace)

#Lets check for reverse polarity and fix 
	finalstream=Stream()
	for trace in st.select(component="Z"):
		dipval = getdip(net,cursta,trace.stats.location,trace.stats.channel,eventtime,sp)
		if debug:
			print 'Here is the dip value:' + str(dipval)
		if dipval == 90.0:
			trace.data = -trace.data
		finalstream += trace
		st.remove(trace)

#Now we rotate the horizontals to E and W
	locations=[]
	for trace in st:
		locations.append(trace.stats.location)
	locations=set(locations)

	for curloc in locations:
		curlochorizontal = st.select(location=curloc)
		curlochorizontal.sort(['channel'])
		if debug:
			print "Here are the number of traces:" + str(len(curlochorizontal)) + " which should be 2"
			print(curlochorizontal)
		azi=getorientation(net,cursta,curloc,curlochorizontal[0].stats.channel,eventtime,sp)
		if debug:
			print "Here is the azimuth" + str(azi)
		curlochorizontal = choptocommon(curlochorizontal)
		finalstream += rotatehorizontal(curlochorizontal,azi)	
			
	if debug:
		print(finalstream)

#We now have rotated data and filtered data so it is time to read in the synthetics and process them
	syns = glob.glob(synfile + '/' + cursta + '.*.LX*.modes.sac*')
	synstream = Stream()
	for cursyn in syns:
		if debug:
			print(cursyn)
		curtrace = read(cursyn)
		curtrace[0].data = curtrace[0].data/(10**9)
		curtrace[0].stats.starttime = curtrace[0].stats.starttime + tshift/2
		if newVer:
			curtrace.taper(max_percentage=0.05, type='cosine')
		else:
			curtrace.taper()
#		pazfake= {'poles': [1-1j], 'zeros': [1-1j], 'gain':1,'sensitivity': 1}
#		curtrace.simulate(paz_remove=pazfake, paz_simulate=pazfake)
		curtrace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
		curtrace.integrate()
		curtrace.integrate()
		if newVer:
			curtrace.taper(max_percentage=0.05, type='cosine')
		else:
			curtrace.taper()
		curtrace.detrend()
		curtrace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
		curtrace[0].stats.channel=(curtrace[0].stats.channel).replace('LH','LX')
		synstream += curtrace
			
	if debug:
		print(synstream)

#We now need to plot everything and save it
#Lets plot the verticals first
	vertcomps = finalstream.select(component="Z")
	vertcomps += synstream.select(component="Z")
	vertcomps = choptocommon(vertcomps)
	if debug:
		print 'Here are the chopped components'
		print(vertcomps)
		
#Set the time series
	tz=numpy.arange(0,vertcomps[0].stats.npts / vertcomps[0].stats.sampling_rate, vertcomps[0].stats.delta)
#Get a legend and plot the vertical
	synplot = figure(1)
	subplot(311)
	titlelegend = vertcomps[0].stats.network + ' ' + vertcomps[0].stats.station + ' '
	stime = str(vertcomps[0].stats.starttime.year) + ' ' + str(vertcomps[0].stats.starttime.julday) + ' ' + \
	str(vertcomps[0].stats.starttime.hour) + ':' + str(vertcomps[0].stats.starttime.minute) + \
	':' + str("{0:.2f}".format(vertcomps[0].stats.starttime.second))
	titlelegend = titlelegend + stime + ' ' 
	lat,lon = getlatlon(cursta,eventtime,sp)
	if debug:
		print "Latitude:" + str(lat)
		print "Longitude:" + str(lon)	
		print "CMT Latitude:" + str(cmtlat)
		print "CMT Longitude:" + str(cmtlon)
	dist= gps2DistAzimuth(float(cmtlat),float(cmtlon),lat,lon)
	bazi ="{0:.1f}".format(dist[2])
	dist ="{0:.1f}".format( 0.0089932 * dist[0] / 1000)
	
	titlelegend = titlelegend + 'Dist:' + str(dist) 
	titlelegend = titlelegend + ' BAzi:' + str(bazi) 
	minper = "{0:.0f}".format(1/usermaxfre)
	maxper = "{0:.0f}".format(1/userminfre)
	titlelegend = titlelegend + ' ' + str(minper) + '-' + str(maxper) + ' s per.'
	title(titlelegend,fontsize=12)
	vertcomps.sort(['location'])
	for comps in vertcomps.select(component="Z"):
		curcolor = getcolor(comps.stats.channel,comps.stats.location)
		plot(tz,(comps.data*(10**3)), curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend(prop={'size':6})


	finalstream += synstream
	finalstream = choptocommon(finalstream)
	finalstream.sort(['location','channel'])
	if debug:
		print "Here is the final stream:"
		print(finalstream)
	subplot(312)
	tne=numpy.arange(0,finalstream[0].stats.npts / finalstream[0].stats.sampling_rate, finalstream[0].stats.delta)
	for comps in finalstream.select(component="N"):
		curcolor = getcolor(comps.stats.channel,comps.stats.location)
		
		plot(tne,(comps.data*(10**3)),curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend(prop={'size':6})
	ylabel('Displacement (mm)')	
	subplot(313)
	for comps in finalstream.select(component="E"):
		curcolor = getcolor(comps.stats.channel,comps.stats.location)
		plot(tne,(comps.data*(10**3)),curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend(prop={'size':6})
	xlabel('Time (s)')
	savefig(os.getcwd() + '/' + resultdir + '/' + cursta + \
	str(vertcomps[0].stats.starttime.year) + str(vertcomps[0].stats.starttime.julday) + \
	str(vertcomps[0].stats.starttime.hour) + str(vertcomps[0].stats.starttime.minute) + '.jpg', format = 'jpeg', dpi=400)

	synplot.clear()

#Time to write some info into the statfile
#Write the network and the station
	#statfile.write(net + "," + cursta + "," + str(dist) + "," + str(bazi))
	writestats(statfile,vertcomps,'Z')
	writestats(statfile,finalstream,'N')
	writestats(statfile,finalstream,'E')
	
	
#Lets get an RMS from the synthetic and the data
	









statfile.close()
