#!/usr/bin/env python

import sys
import os
import glob
import numpy as np
import matplotlib
import math
import warnings
import argparse

import matplotlib.pyplot as plt
from obspy import read, Stream
from obspy.core import UTCDateTime
from obspy.io.xseed import Parser

from time import gmtime, strftime
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.cross_correlation import xcorr
from obspy.core.util.version import read_release_version



def getorientation(tr, sp):
    sta, net, chan, loc = getsncl(tr)
    evetime = tr.stats.starttime
    # A function to get the orientation of a station at a specific time
    for cursta in sp.stations:
        # As we scan through blockettes we need to find blockettes 50 and 52
        for blkt in cursta:
            if blkt.id == 50:
                # Pull the station info for blockette 50
                stacall = blkt.station_call_letters.strip()
            if stacall == sta:
                if blkt.id == 52 and blkt.location_identifier == loc and blkt.channel_identifier == chan:
                    if type(blkt.end_date) is str:
                        curdoy = strftime("%j", gmtime())
                        curyear = strftime("%Y", gmtime())
                        curtime = UTCDateTime(curyear + "-" +
                                              curdoy + "T00:00:00.0")
                        if blkt.start_date <= evetime:
                            azimuth = blkt.azimuth
                    elif blkt.start_date <= evetime and blkt.end_date >= evetime:
                        azimuth = blkt.azimuth
    return azimuth


def getdip(tr, sp):
    sta, net, chan, loc = getsncl(tr)
    evetime = tr.stats.starttime
    # A function to get the dip of a station at a specific time
    dip = 0.
    for cursta in sp.stations:
        # As we scan through blockettes we need to find blockettes 50 and 52
        for blkt in cursta:
            if blkt.id == 50:
                # Pull the station info for blockette 50
                stacall = blkt.station_call_letters.strip()
            if stacall == sta:
                if blkt.id == 52 and blkt.location_identifier == loc and blkt.channel_identifier == chan:
                    if type(blkt.end_date) is str:
                        curdoy = strftime("%j", gmtime())
                        curyear = strftime("%Y", gmtime())
                        curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
                        if blkt.start_date <= evetime:
                            dip = blkt.dip
                    elif blkt.start_date <= evetime and blkt.end_date >= evetime:
                        dip = blkt.dip
                    
    return dip


def rotatehorizontal(stream, angle1, angle2):
    # Switch to E and N
    debugRot = False
    if stream[0].stats.channel in set(['LHE', 'LHN', 'BHE', 'BHN']):
        stream.sort(['channel'], reverse=True)
        angle1, angle2 = angle2, angle1
    if debugRot:
        print(stream)
        print 'Angle1: ' + str(angle1) + ' Angle2: ' + str(angle2)
    theta_r1 = math.radians(angle1)
    theta_r2 = math.radians(angle2)
    swapSecond = False
    if (angle2 >= 180. and angle2 <= 360.) or angle2 == 0.:
        swapSecond = True 
    # if the components are swaped swap the matrix
    if theta_r1 > theta_r2 and swapSecond:
        if debugRot:
            print 'Swap the components: ' + str((360. - angle1) - angle2)
        stream.sort(['channel'], reverse=True)
        theta_r1, theta_r2 = theta_r2, theta_r1
        print(stream)
    # create new trace objects with same info as previous
    rotatedN = stream[0].copy()
    rotatedE = stream[1].copy()
    # assign rotated data
    rotatedN.data = stream[0].data*math.cos(-theta_r1) +\
        stream[1].data*math.sin(-theta_r1)
    rotatedE.data = -stream[1].data*math.cos(-theta_r2-math.pi/2.) +\
        stream[0].data*math.sin(-theta_r2-math.pi/2.)
    rotatedN.stats.channel = 'LHN'
    rotatedE.stats.channel = 'LHE'
    # return new streams object with rotated traces
    streamsR = Stream(traces=[rotatedN, rotatedE])
    return streamsR


def choptocommon(stream):
    # A function to chop the data to a common time window
    stime = max([tr.stats.starttime for tr in st])
    etime = min([tr.stats.endtime for tr in st])
    stream.trim(starttime=stime, endtime=etime)
    return stream

def getstalist(sp, etime, net):
    # A function to get a station list
    stations = []
    for cursta in sp.stations:
        # As we scan through blockettes we need to find blockettes 50 
        for blkt in cursta:
            if blkt.id == 50:
                # Pull the station info for blockette 50
                stacall = blkt.station_call_letters.strip()
                if debug:
                    print "Here is a station in the dataless: " + stacall
                if type(blkt.end_effective_date) is str:
                    curdoy = strftime("%j", gmtime())
                    curyear = strftime("%Y", gmtime())
                    curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
                    
                    if blkt.start_effective_date <= etime:
                        stations.append(blkt.station_call_letters.strip())
                elif blkt.start_effective_date <= etime and blkt.end_effective_date >= etime:
                    stations.append(blkt.station_call_letters.strip())    
    return stations


def readcmt(cmt, debug=False):
    # Now we can continue like there is no difference between 
    # Princeton and our Synthetics
    # Lets get the event time from the cmt
    cmtline1 = ' '.join(cmt[0].split())
    cmtlat = cmt[4].replace('latitude:', '').strip()
    cmtlon = cmt[5].replace('longitude:', '').strip()
    tshift = float(cmt[2].replace('time shift:', '').strip())
    hdur = float(cmt[3].replace('half duration:', '').strip())
    if debug:
        print(cmtline1)
    cmtline1 = cmtline1.split()
    if debug:
        print(cmtline1[1] + ' ' + cmtline1[2] + ' ' +\
              cmtline1[3] + ' ' + cmtline1[4] + ' ' + \
              cmtline1[5] + ' ' + cmtline1[6])
    eventtime = UTCDateTime(int(cmtline1[1]), int(cmtline1[2]),
                            int(cmtline1[3]), int(cmtline1[4]),
                            int(cmtline1[5]), float(cmtline1[6]))
    if debug:
        print 'Year:' + str(eventtime.year)
        print 'Day:' + str(eventtime.julday)
        print 'Hour:' + str(eventtime.hour)
        print 'Minute:' + str(eventtime.minute)
    return cmtlat, cmtlon, eventtime, tshift, hdur


def getdata(net, sta, eventtime, lents, debug=False):
    # This function goes to both archives and gets the data
    stime = eventtime - 3000.
    etime = eventtime + 3000. + lents

    # Current default is LH add BH later
    chan = 'LH'

    if (net in set(['GT'])) or (sta == 'KBL'):
        chan = 'BH'

    # Grab the data locations
    datalocs = ['/tr1/telemetry_days/', '/xs*/seed/']

    files = []
    # Grab the files and deal with edge cases
    for dataloc in datalocs:
        for year in range(stime.year, etime.year+1):
            for day in range(stime.julday, etime.julday+1):
                string = dataloc + net + '_' + sta + '/' + \
                        str(year) + '/' + str(stime.year) + '_' + \
                        str(day).zfill(3) + '*/*' + chan + '*.seed'
                files += glob.glob(string)
    if debug:
        print(files)
    st = Stream()
    for curfile in files:
        try:
            st += read(curfile, starttime=stime, endtime=etime)
        except:
            if debug:
                print 'Unable to get data ' + curfile
    st.merge(fill_value='latest')
    if debug:
        print 'We have data'
    return st
    
def fixdip(st, eventtime, sp, debug=False):
    for tr in st.select(component="Z"):
        dipval = getdip(tr, sp)
        if debug:
            print 'Here is the dip value:' + str(dipval)
        if dipval == 90.0:
            tr.data = -tr.data
    return st


def getcolor(chan, loc):
    if chan in set(['LXN', 'LXE', 'LXZ']):
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


def writestats(statfile, streamin, comp):
    try:
        syncomp = "LX" + comp    
        datacomp = "LH" + comp
        syn = streamin.select(channel = syncomp)
        for tr in streamin.select(channel = datacomp):    
            resi = "{0:.2f}".format(np.sum(tr.data*syn[0].data)/np.sum(np.square(syn[0].data)))
            lag, corr = xcorr(tr,syn[0],500)
            corr = "{0:.2f}".format(corr)
            statfile.write(tr.stats.network + "," + tr.stats.station)
            statfile.write("," + tr.stats.location + "," + tr.stats.channel + "," +  str(resi))
            statfile.write("," + str(lag) + "," + str(corr) + ", ")
            statfile.write(str(tr.stats.starttime.month) + "/" + str(tr.stats.starttime.day) + \
                "/" + str(tr.stats.starttime.year) + " " + str(tr.stats.starttime.hour) + ":" + \
                str(tr.stats.starttime.minute) + ":" + str(tr.stats.starttime.second) + "\n")
    except:    
        if debug:
            print 'No residual for' + tr.stats.station + ' ' + 'LH' + comp    
    return


def getargs():
    # Grab command line arguments to run synthetics
    parser = argparse.ArgumentParser(description = "Program to compare long-period event synthetics to data")

    parser.add_argument('-n', type=str, action="store",
                        dest = "network", required=True,
                        help="Network name Example: IU")

    parser.add_argument('-resDir',type=str, action = "store",
                        dest="resDir", required=True,
                        help="Result directory name Example: blah")

    parser.add_argument('-syn', type=str, action="store",
                        dest="syn", required=True, nargs="+",
                        help="Synthetics directory location Example: " + \
                        "/SYNTHETICS/2014/C201401*")

    parser.add_argument('-sta', type=str, action="store",
                        dest="sta", required = False,
                        help="Stations to use Example with a comma (,) separator : TUC,ANMO")

    parser.add_argument('-tslen', type=int, action="store",
                        dest="lents", required=False, default=4000.,
                        help="Length of time series in seconds Example:  4000, default is 4000 s")

    parser.add_argument('-debug', action="store_true", dest="debug",
        default=False, help="Run in debug mode")

    parser.add_argument('-filter', action="store", nargs=3, dest="filter", required=False,
        default = [ 100, 200, 4], 
        help="Filter parameters using minimum period maximum period and number of corners Example: 100 200 4, " + \
            "default is 100 400 2")

    parserval=parser.parse_args()
    return parserval

    
def procStream(st, sp, eventtime, tshift, freqmin, freqmax, corners, lents, debug=False):
    # Lets go through each trace in the stream and deconvolve and filter
    if len(st.select(channel="LX*")) > 0:
        synthetic = True
    else:
        synthetic = False
    if debug:
        print(st)
        print('Synthetics: ' + str(synthetic))
    for tr in st:
        try:
            # Here we get the response and remove it
            # Here is where I am mucking around
            if synthetic:
                tr.data /= (10**9)
                tr.stats.starttime += float(tshift)/2.
            else:
                paz=sp.get_paz(tr.id, eventtime)
                if debug:
                    print 'Here is the paz'
                    print(paz)        
        except:
            print('Problem with response')
            st.remove(tr)
        try:    
            tr.taper(max_percentage=0.05, type='cosine')
            if not synthetic:
                tr.simulate(paz_remove=paz)
            # Here we filter
            tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=corners)
            tr.integrate()
            if synthetic:
                tr.integrate()
            tr.taper(max_percentage=0.05, type='cosine')
            tr.trim(starttime=eventtime + float(tshift)/2.,endtime=(eventtime+lents + float(tshift)/2.))
            tr.detrend()
            tr.filter("bandpass", freqmin=freqmin,freqmax=freqmax, corners=corners)
            if synthetic:
                tr.stats.channel = (tr.stats.channel).replace('LH','LX')
        except:
            print('Problem filtering')
            st.remove(tr)
    
    return st
    
def pltStream(stream, pltHandle, component, cmtlat=None, cmtlon=None,
              minfre=None, maxfre=None, resDir=None, debug = False):
    if component == 'Z':
        plt.subplot(3,1,1)
        title = stream[0].stats.network + ' ' + \
                stream[0].stats.station
        starttime = stream[0].stats.starttime
        stime = str(starttime.year) + ' ' + str(starttime.julday) + \
                ' ' + str(starttime.hour) + ':' + \
                str(starttime.minute) + ':' + \
                str("{0:.2f}".format(starttime.second))
        title += stime + ' '
        coord = sp.get_coordinates(stream[0].id, starttime)
        lat = coord['latitude']
        lon = coord['longitude']
        if debug:
            print "Latitude:" + str(lat)
            print "Longitude:" + str(lon)
            print "CMT Latitude:" + str(cmtlat)
            print "CMT Longitude:" + str(cmtlon)
        dist = gps2dist_azimuth(float(cmtlat), float(cmtlon), lat, lon)
        bazi = "{0:.1f}".format(dist[2])
        dist = "{0:.1f}".format(0.0089932*dist[0]/1000.)
        title += 'Dist:' + str(dist)
        title += ' BAzi:' + str(bazi)
        # Add period band
        title += str("{0:.0f}".format(1/maxfre)) + '-' + \
                 str("{0:.0f}".format(1/minfre)) + ' s per.'
        # The title is now finished
        plt.title(title, fontsize=12)
    elif component == 'N':
        plt.subplot(3,1,2)
        plt.ylabel('Displacement (mm)')
    elif component == 'E':
        plt.subplot(3,1,3)
        plt.xlabel('Time (s)')
    else:
        print('Unable to plot: ' + str(component))
        return
    t = np.arange(0, stream[0].stats.npts)/stream[0].stats.sampling_rate
    stream.sort(['location'])
    for tr in stream.select(component=component):
        color = getcolor(tr.stats.channel, tr.stats.location)
        plt.plot(t, tr.data*(10**3), color, label = tr.stats.location +
                 ' ' + tr.stats.channel)
    plt.legend(prop={'size': 6}, loc=2)
    plt.xlim((min(t), max(t)))
    if component == 'E':
        starttime = stream[0].stats.starttime
        plt.savefig(os.getcwd() + '/' + resDir + '/' + 
                    stream[0].stats.network + stream[0].stats.station +
                    str(starttime.year) + str(starttime.julday) +
                    str(starttime.hour) + str(starttime.minute) +
                    '.jpg', format='jpeg', dpi=400)
        synplot.clear()
    return

def getsncl(tr):
    # Return the snlc
    nslc = (tr.id).split('.')
    return nslc[1], nslc[0], nslc[3], nslc[2]
                


# Start of the main part of the program
if __name__ == "__main__":
    
    datalessloc = '/APPS/metadata/SEED/'

    # Lets get the parser arguments
    parserval = getargs()

    # Debug flag
    debug = parserval.debug

    # Length of time series
    lents = parserval.lents

    # Grab the filter parameters
    userminfre = 1.0/float(parserval.filter[1])
    usermaxfre = 1.0/float(parserval.filter[0])
    filtercornerpoles = int(parserval.filter[2])

    # Lets read in the dataless
    net = parserval.network
    try:
        sp = Parser(datalessloc + net + ".dataless")
    except:
        print "Can not read the dataless."
        exit(0)

    # Run through each of the event CMTs
    for synfile in parserval.syn:
     
        if synfile[-1] == '/':
            synfile = synfile[:-1]

        # Read in the CMT solution         
        if debug:
            print "We are using local synthetics"
        if not os.path.isfile(synfile + '/CMTSOLUTION'):
            print "No CMT found"
            exit(0)
        cmt = tuple(open(synfile + '/CMTSOLUTION'))
        cmtlat, cmtlon, eventtime, tshift, hdur = readcmt(cmt)


        # Lets make a local results directory
        resultdir = parserval.resDir
        if resultdir[-1] == '/':
            resultdir = resultdir[:-1]

        if not os.path.exists(os.getcwd() + '/' + resultdir):
            os.mkdir(os.getcwd() + '/' + resultdir)
        evename = synfile.split("/")
        evename = evename[len(evename)-1]

        statfile = open(os.getcwd() + '/' + resultdir + '/Results' + evename + net + '.csv' ,'w')
        statfile.write('net,sta,loc,chan,scalefac,lag,corr,time\n')

            
        if parserval.sta:
            if debug: 
                print "We are using a manual station list"
            stations = parserval.sta.split(",")
        else:
            stations = getstalist(sp, eventtime, net)

        if debug:
            print "Here are the stations we found"    
            for sta in stations:
                print "Here is a station:" + sta


        # Lets start by using a station list and then move to a different approach
        for sta in stations:
            #try:
            if True:
                sta = sta.strip()
                st = getdata(net, sta, eventtime, lents)
                if len(st) == 0:
                    pass
                    print 'Bad data for: ' + sta
                # Decimate any high sample rate traces to 1
                for tr in st:
                    if tr.stats.sampling_rate > 20:
                            tr.detrend('linear')
                            tr.taper(max_percentage=0.05, type='cosine')
                            tr.decimate(int(tr.stats.sampling_rate/4.))
                            tr.decimate(4)

                # Lets go through each trace in the stream and deconvolve and filter
                st = procStream(st, sp, eventtime, tshift, userminfre, usermaxfre, filtercornerpoles, lents)
                # Lets check for reverse polarity and fix

                st = fixdip(st, eventtime, sp)
                stF = st.select(component="Z")
                for tr in st.select(component="Z"):
                    st.remove(tr)
                # Now we rotate the horizontals to E and W
                locations = set([tr.stats.location for tr in st])

                for curloc in locations:
                    curlochorizontal = st.select(location=curloc)
                    curlochorizontal.sort(['channel'])
                    if debug:
                        print "Here are the number of traces:" + str(len(curlochorizontal)) + \
                            " which should be 2"
                        print(curlochorizontal)
                    azi1=getorientation(curlochorizontal[0], sp)
                    azi2=getorientation(curlochorizontal[1], sp)   
                    if debug:
                        print "Here is the azimuth" + str(azi1) + ' ' + str(azi2)
                    curlochorizontal = choptocommon(curlochorizontal)
                    stF += rotatehorizontal(curlochorizontal, azi1, azi2)    
            
                for tr in stF.select(channel="BH*"):
                    tr.stats.channel = tr.stats.channel.replace('B','L')

                files = glob.glob(synfile + '/' + sta + '.*.LX*.modes.sac')
                synstream = Stream()
                for curfile in files:
                   synstream += read(curfile)
                for tr in synstream:
                    tr.stats.channel = (tr.stats.channel).replace('LH','LX')
                    tr.stats.network = net  
                synstream = procStream(synstream, sp, eventtime, tshift, userminfre, usermaxfre, filtercornerpoles, lents)
                
                stF += synstream

                synplot =  plt.figure(1)
                for comp in ["Z", "N", "E"]:    
                    pltStream(stF, synplot, comp,
                              cmtlat=cmtlat, cmtlon=cmtlon,
                              minfre=userminfre, maxfre=usermaxfre,
                              resDir=resultdir)
                    # Time to write some info into the statfile
                    writestats(statfile, stF, comp)
        statfile.close()
