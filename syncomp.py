#!/usr/bin/env python

import sys
import os
import glob
import numpy
import matplotlib
import math
from matplotlib.pyplot import (figure,axes,plot,xlabel,ylabel,title,subplot,legend,savefig,show)
from obspy import read, Stream
from obspy.core import UTCDateTime
from obspy.xseed import Parser
from time import gmtime, strftime
from obspy.core.util.geodetics import gps2DistAzimuth


debug=False
datalessloc = '/APPS/metadata/SEED/'
userminfre = .005
usermaxfre = .01
lents = 8000


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
#Now we can continue like there is no difference between Princeton and our Synthetics
#Lets get the event time from the cmt
	cmtline1 = ' '.join(cmt[0].split())
	cmtlat = cmt[4].replace('latitude:','').strip()
	cmtlon = cmt[5].replace('longitude:','').strip()
	if debug:
		print cmtline1
	cmtline1 = cmtline1.split()
	if debug:
		print cmtline1[1] + ' ' + cmtline1[2] + ' ' + cmtline1[3] + ' ' + cmtline1[4] + ' ' + cmtline1[5] + ' ' + cmtline1[6]
	eventtime = UTCDateTime(int(cmtline1[1]),int(cmtline1[2]),int(cmtline1[3]),int(cmtline1[4]),int(cmtline1[5]),float(cmtline1[6]))
	if debug:
		print 'Year:' + str(eventtime.year)
		print 'Day:' + str(eventtime.julday)
		print 'Hour:' + str(eventtime.hour)
		print 'Minute:' + str(eventtime.minute)
	return cmtlat, cmtlon, eventtime


def getdata(net,sta,eventtime,lents):
	if net in set(['IW','NE','US']):	
		dataloc = 'xs1'	
	elif net in set(['IU','IC','CU']):
		dataloc = 'xs0'	
	
	st = read('/' + dataloc + '/seed/' + net + '_' + sta + '/' + str(eventtime.year) + \
	'/' + str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '_' + net + '_' + sta + '/*LH*.seed', \
	starttime=eventtime-3000,endtime=(eventtime+lents+3000))
	st += read('/' + dataloc + '/seed/' + net + '_' + sta + '/' + str(eventtime.year) + \
	'/' + str(eventtime.year) + '_' + str(eventtime.julday + 1).zfill(3) + '_' + net + '_' + sta + '/*LH*.seed', \
	starttime=eventtime-3000,endtime=(eventtime+lents+3000))
	st += read('/' + dataloc + '/seed/' + net + '_' + sta + '/' + str(eventtime.year) + \
	'/' + str(eventtime.year) + '_' + str(eventtime.julday - 1).zfill(3) + '_' + net + '_' + sta + '/*LH*.seed', \
	starttime=eventtime-3000,endtime=(eventtime+lents+3000))
	st.merge(fill_value='latest')
	if debug:
		print 'We have data'

	return st




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
cmtlat, cmtlon, eventtime = readcmt(cmt)

#Lets make a local results directory
resultdir = sys.argv[2]
if not os.path.exists(os.getcwd() + '/' + resultdir):
	os.mkdir(os.getcwd() + '/' + resultdir)

curnet = sys.argv[3]

#Lets read in the dataless
try:
	sp = Parser(datalessloc + curnet + ".dataless")
except:
	print "Can not read the dataless."
	exit(0)

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
		st = getdata(net,cursta,eventtime,lents)
	except:
		print('No data for ' + net + ' ' + cursta)
		continue
		
#Lets go through each trace in the stream and deconvolve and filter
	for trace in st:
#Here we get the response and remove it
		try:
			paz=sp.getPAZ(net + '.' + cursta + '.' + trace.stats.location + '.' + trace.stats.channel,datetime=eventtime)
			trace.integrate()
			trace.taper(type='cosine')				
			trace.simulate(paz_remove=paz)
#Here we filter
			trace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=4)
			trace.taper(type='cosine')
			trace.trim(starttime=eventtime,endtime=(eventtime+lents))
		except:
			print('Can not find the response')
			st.remove(trace)
#Now we rotate the horizontals to E and W
	horizontalstream = st	
	finalstream=Stream()
	for trace in horizontalstream.select(component="Z"):
		dipval = getdip(net,cursta,trace.stats.location,trace.stats.channel,eventtime,sp)
		if debug:
			print 'Here is the dip value:' + str(dipval)
		if dipval == 90:
			trace.data = -trace.data
		finalstream += trace
		horizontalstream.remove(trace)

	if debug:
		print 'Here are the horizontal traces:'
		print(horizontalstream)

	locations=[]
	for trace in horizontalstream:
		locations.append(trace.stats.location)
	locations=set(locations)

	for curloc in locations:
		curlochorizontal = horizontalstream.select(location=curloc)
		if debug:
			print "Here are the number of traces:" + str(len(curlochorizontal)) + " which should be 2"
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
		curtrace.integrate()
		curtrace.integrate()
		curtrace[0].data = curtrace[0].data/(10**9)
		curtrace.taper(type='cosine')
		pazfake= {'poles': [1-1j], 'zeros': [1-1j], 'gain':1,'sensitivity': 1}
		curtrace.simulate(paz_remove=pazfake, paz_simulate=pazfake)
		curtrace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=4)
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
	dist ="{0:.2f}".format( 0.0089932 * dist[0] / 1000)
	titlelegend = titlelegend + 'Distance:' + str(dist) + ' degrees'
	title(titlelegend)
	vertcomps.sort(['location'])
	for comps in vertcomps:
		if comps.stats.channel == 'LXZ':
			curcolor = 'k'
		elif comps.stats.location == '00':
			curcolor = 'g'
		elif comps.stats.location == '10':
			curcolor = 'r'
		else:
			curcolor = 'c'
		plot(tz,(comps.data*(10**3)), curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend()


	finalstream += synstream.select(component="N")
	finalstream += synstream.select(component="E")
	finalstream = choptocommon(finalstream)
	finalstream.sort(['location','channel'])
	if debug:
		print "Here is the final stream:"
		print(finalstream)
	subplot(312)
	tne=numpy.arange(0,finalstream[0].stats.npts / finalstream[0].stats.sampling_rate, finalstream[0].stats.delta)
	for comps in finalstream.select(component="N"):
		if comps.stats.channel == 'LXN':
			curcolor = 'k'
		elif comps.stats.location == '00':
			curcolor = 'g'
		elif comps.stats.location == '10':
			curcolor = 'r'
		else:
			curcolor = 'c'
		plot(tne,(comps.data*(10**3)),curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend()
	ylabel('Displacement (mm)')	
	subplot(313)
	for comps in finalstream.select(component="E"):
		if comps.stats.channel == 'LXE':
			curcolor = 'k'
		elif comps.stats.location == '00':
			curcolor = 'g'
		elif comps.stats.location == '10':
			curcolor = 'r'
		else:
			curcolor = 'c'
		plot(tne,(comps.data*(10**3)),curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend()
	xlabel('Time (s)')
	savefig(os.getcwd() + '/' + resultdir + '/' + cursta + \
	str(vertcomps[0].stats.starttime.year) + str(vertcomps[0].stats.starttime.julday) + \
	str(vertcomps[0].stats.starttime.hour) + str(vertcomps[0].stats.starttime.minute) + '.jpg', format = 'jpeg')

	synplot.clear()


