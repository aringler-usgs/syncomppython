#!/usr/bin/env python

import sys
import os
import glob
import numpy
import matplotlib
import math
import warnings
import argparse
from matplotlib.pyplot import (figure,axes,plot,xlabel,ylabel,title,subplot,legend,savefig,show)
from obspy import read, Stream
from obspy.core import UTCDateTime
from obspy.xseed import Parser
from time import gmtime, strftime
from obspy.core.util.geodetics import gps2DistAzimuth
from obspy.signal.cross_correlation import xcorr
from obspy.core.util.version  import read_release_version


datalessloc = '/APPS/metadata/SEED/'


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
	if stream[0].stats.channel in set(['LHE','LHN','BHE','BHN']):
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
					print "Here is a station in the dataless: " + stacall
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
		print cmtline1[1] + ' ' + cmtline1[2] + ' ' + cmtline1[3] + ' ' + cmtline1[4] + ' ' + \
			cmtline1[5] + ' ' + cmtline1[6]
	eventtime = UTCDateTime(int(cmtline1[1]),int(cmtline1[2]),int(cmtline1[3]),int(cmtline1[4]),\
		int(cmtline1[5]),float(cmtline1[6]))
	if debugreadcmt:
		print 'Year:' + str(eventtime.year)
		print 'Day:' + str(eventtime.julday)
		print 'Hour:' + str(eventtime.hour)
		print 'Minute:' + str(eventtime.minute)
	return cmtlat, cmtlon, eventtime, tshift,hdur

def getdata(net,sta,eventtime,lents,dataloc):
#This function goes to one of the archives and gets the data
	debuggetdata = False
	preeventday = eventtime - 24*60*60
	posteventday = eventtime + 24*60*60
	prepostwin= 3000
	chanType = 'LH'
	if (net in set(['GT'])) or (cursta == 'KBL'):
		chanType = 'BH'
#If II get off of /tr1 else get the data from /xs0 or /xs1
	if net == 'II':
		dataprefix1 = '/tr1/telemetry_days/'
		dataperfix = [dataprefix1]
	else:
		if net in set(['IW','NE','US']):	
			dataprefix = 'xs1'	
		else:
			dataprefix = 'xs0'
		dataprefix1 = '/' + dataprefix + '/seed/'
		dataprefix2 = '/tr1/telemetry_days/'
		if dataloc:
			dataprefix = [dataprefix1]
		else:
			dataprefix = [dataprefix1, dataprefix2]
	if debuggetdata:
		print 'Here is the dataprefix:' + dataprefix1
	st = Stream()
	for dataprefixs in dataprefix:
		try:
			st += read(dataprefixs + net + '_' + sta + '/' + str(eventtime.year) + \
				'/' + str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '*/*' + chanType + '*.seed', \
				starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
			st += read(dataprefixs + net + '_' + sta + '/' + str(posteventday.year) + \
				'/' + str(posteventday.year) + '_' + str(posteventday.julday).zfill(3) + '*/*' + chanType + '*.seed', \
				starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
			st += read(dataprefixs + net + '_' + sta + '/' + str(preeventday.year) + \
				'/' + str(preeventday.year) + '_' + str(preeventday.julday).zfill(3) + '*/*' + chanType + '*.seed', \
				starttime=eventtime-prepostwin,endtime=(eventtime+lents+prepostwin))
		except:
			if debuggetdata:
				print 'Unable to get data'

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
				elif blockette.stage_sequence_number == 3:
					if not 'digitizer_gain' in data.keys():
						data['digitizer_gain'] = blockette.sensitivity_gain
					else:
						data['digitizer_gain'] *= blockette.sensitivity_gain
			elif blockette.id == 53 and channel_flag and station_flag:
				if not 'gain' in data.keys():
					data['gain'] = blockette.A0_normalization_factor
				else:
					data['gain'] *= blockette.A0_normalization_factor
				if debuggetPAZ2:
					print 'Here is the gain: ' + str(blockette.A0_normalization_factor)
				if not 'poles' in data.keys():
					data['poles'] = []
				if not blockette.transfer_function_types in set(['A','B']):
					msg = 'Only supporting Laplace transform response ' + \
					'type. Skipping other response information.'
					warnings.warn(msg, UserWarning)
					continue
				
				if blockette.transfer_function_types == 'B':
					schan = 1
				else:
					schan = 1
				if debuggetPAZ2:
					print 'Here are the number of poles:' + str(blockette.number_of_complex_poles)
					print 'Here are the number of zeros:' + str(blockette.number_of_complex_zeros)
				for i in range(blockette.number_of_complex_poles):
					p = complex(schan*blockette.real_pole[i], schan*blockette.imaginary_pole[i])
					data['poles'].append(p)
				if not 'zeros' in data.keys():
					data['zeros'] = []
				for i in range(blockette.number_of_complex_zeros):
					z = complex(schan*blockette.real_zero[i], schan*blockette.imaginary_zero[i])
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
	elif loc == '30':
		color = '0.75'
	elif loc == '40':
		color = 'y'
	elif loc == '50':
		color = 'm'
	else:
		color = 'b'
	return color


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


def getargs():

	parser = argparse.ArgumentParser(description = "Program to compare long-period event synthetics to data")

	parser.add_argument('-n', type = str, action = "store", \
		dest = "network", required = True, help = "Network name Example: IU")

	parser.add_argument('-resDir',type = str, action = "store", \
		dest = "resDir", required = True, help = "Result directory name Example: blah")

	parser.add_argument('-syn', type = str, action = "store", \
		dest = "syn", required = True, nargs = "+", help = "Synthetics directory location Example: " + \
		"/SYNTHETICS/2014/C201401*")

	parser.add_argument('-sta', type = str, action = "store", \
		dest = "sta", required = False, help = "Stations to use Example with a comma (,) separator : TUC,ANMO")

	parser.add_argument('-tslen', type = int, action ="store", \
		dest = "lents", required = False, help = "Length of time series in seconds Example:  4000, default is 4000 s")

	parser.add_argument('-debug', action = "store_true", dest = "debug", \
		default = False, help = "Run in debug mode")

	parser.add_argument('-dataloc', action = "store_true", dest = "dataloc", \
		default = False, help = "Use /xs0 data location, otherwise use /tr1 also")

	parser.add_argument('-filter', action = "store", nargs = 3, dest = "filter", required = False, \
		help = "Filter parameters using minimum period maximum period and number of corners Example: 100 200 4, " + \
			"default is 100 400 2")

	parserval = parser.parse_args()
	return parserval




#Start of the main part of the program
if __name__ == "__main__":

	#Lets get the parser arguments
	parserval = getargs()

	debug = parserval.debug

	if parserval.sta:
		manstalist = True
		if debug: 
			print "We are using a manual station list"
		stalist = parserval.sta.split(",")
		stations = []
		for sta in stalist:
			stations.append(parserval.network + " " + sta)
		if debug:
			print(stations) 
	else:
		manstalist = False

	if parserval.lents:
		lents = parserval.lents
	else:
		lents = 4000

	if parserval.filter:
		userminfre = 1.0/float(parserval.filter[1])
		usermaxfre = 1.0/float(parserval.filter[0])
		filtercornerpoles = int(parserval.filter[2])
	else:
		userminfre = 0.0025
		usermaxfre = 0.01
		filtercornerpoles = 2

	if parserval.dataloc:
		dataloc = True
	else:
		dataloc = False


	#Lets read in the dataless
	curnet = parserval.network
	try:
		sp = Parser(datalessloc + curnet + ".dataless")
	except:
		print "Can not read the dataless."
		exit(0)

	for synfile in parserval.syn:
	 

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
		resultdir = parserval.resDir
		if resultdir[-1] == '/':
			resultdir = resultdir[:-1]

		if not os.path.exists(os.getcwd() + '/' + resultdir):
			os.mkdir(os.getcwd() + '/' + resultdir)
		evename = synfile.split("/")
		evename = evename[len(evename)-1]

		

		statfile = open(os.getcwd() + '/' + resultdir + '/Results' + evename + curnet + '.csv' ,'w')
		statfile.write('net,sta,loc,chan,scalefac,lag,corr\n')

		
		if not manstalist:
			stations = getstalist(sp,eventtime,curnet)

		if debug:
			print "Here are the stations we found"	
			for sta in stations:
				print "Here is a station:" + sta


		#Lets start by using a station list and then move to a different approach
		for sta in stations:
			try:
				cursta = sta.strip()
				if debug:
					print 'Current station:' + cursta
				cursta = sta.split()
				net = cursta[0]
				cursta = cursta[1]
		
				st = getdata(net,cursta,eventtime,lents,dataloc)
				if len(st) == 0:
					pass
					print 'Bad data for: ' + sta
				if (net in set(['GT'])) or (cursta == 'KBL'):
					for tr in st:
						if tr.stats.sampling_rate > 20:
							tr.detrend('linear')
							tr.taper(max_percentage=0.05, type='cosine')
							tr.decimate(int(tr.stats.sampling_rate/4.))
							tr.decimate(4)

				#Lets go through each trace in the stream and deconvolve and filter
				for trace in st:
					#Here we get the response and remove it
					#Here is where I am mucking around		
					paz=getPAZ2(sp,net,cursta,trace.stats.location,trace.stats.channel,eventtime)
					if debug:
						print 'Here is the paz'
						print(paz)		
					try:
			
						trace.taper(max_percentage=0.05, type='cosine')
						trace.simulate(paz_remove=paz)
						#Here we filter
						trace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
						trace.integrate()
						trace.taper(max_percentage=0.05, type='cosine')

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
						print "Here are the number of traces:" + str(len(curlochorizontal)) + \
							" which should be 2"
						print(curlochorizontal)
					azi=getorientation(net,cursta,curloc,curlochorizontal[0].stats.channel,eventtime,sp)
					if debug:
						print "Here is the azimuth" + str(azi)
					curlochorizontal = choptocommon(curlochorizontal)
					finalstream += rotatehorizontal(curlochorizontal,azi)	
			
				if debug:
					print(finalstream)
				if (net in set(['GT'])) or (cursta == 'KBL'):
					for tr in finalstream:
						tr.stats.channel = tr.stats.channel.replace('B','L')
			#We now have rotated data and filtered data so it is time to read in the synthetics and process them
				syns = glob.glob(synfile + '/' + cursta + '.*.LX*.modes.sac')
				synstream = Stream()
				for cursyn in syns:
					if debug:
						print(cursyn)
					curtrace = read(cursyn)
					curtrace[0].data = curtrace[0].data/(10**9)
					curtrace[0].stats.starttime = curtrace[0].stats.starttime + tshift/2

					curtrace.taper(max_percentage=0.05, type='cosine')

			#		pazfake= {'poles': [1-1j], 'zeros': [1-1j], 'gain':1,'sensitivity': 1}
			#		curtrace.simulate(paz_remove=pazfake, paz_simulate=pazfake)
					curtrace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
					curtrace.integrate()
					curtrace.integrate()
		
					curtrace.taper(max_percentage=0.05, type='cosine')

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
	
			except:
				print 'Bad data for station: ' + sta




		statfile.close()
