#!/usr/bin/env python

#########################################################################
# syncomp.py
#
# created by Adam Ringler
#   updated some of the obspy import statements
#   added documentation for user clarity
# 
#
# - Adam has a remark about using something other than a station list
#########################################################################

""" 
    syncomp.py 
        Plots up synthetic seismograms calculated from Mineos with data
        as a tool to find instrument problems not easily seen by eye.
        
        to run:
        syncomp.py -syn /SYNTHETICS/2017/C201701032152A -n CU -resDir 2017_003
        where
        -syn /SYNTHETICS/2017/C201701032152A is the directory containing the 
            synthetics and CMT information
        -resDir 2017_003 is the directory to place the plots in
        -n CU is the network you wish to plot 
        -sta GRGR is for a single station plot (do not specify station
            if you want to plot up the entire network) 
"""

import sys
import os
import glob
import numpy as np
import matplotlib
import math
import warnings
import argparse

import matplotlib.pyplot as plt
from obspy import read, Stream, read_events
from obspy.core import UTCDateTime
from obspy.io.xseed import Parser
from obspy.clients.fdsn import Client

from time import gmtime, strftime
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.cross_correlation import xcorr
from obspy.core.util.version import read_release_version
from shutil import copyfile

from scipy import signal




def choptocommon(st):
    """ A function to chop the data to a common time window. """
    stime = max([tr.stats.starttime for tr in st])
    etime = min([tr.stats.endtime for tr in st])
    st.trim(starttime=stime, endtime=etime)
    if debug:
        print('starttime: '+str(stime)+' endtime: '+str(etime))
    return st


def readcmt(cmt, debug=False):
    """ read the CMT information contained in event object """
    hdur = cat[0].focal_mechanisms[0].moment_tensor.source_time_function.duration 
    cmtlat = cat[0].origins[0].latitude
    cmtlon = cat[0].origins[0].longitude
    eventtime = cat[0].origins[0].time
    # The time shift is included in the event time for the cmt
    # This should eventually be removed
    tshift = cat[0].origins[0].time - cat[0].origins[1].time
    return cmtlat, cmtlon, eventtime, tshift, hdur


def getdata(net, sta, eventtime, lents, debug=False):
    """
        This function goes to both archives and gets the data.
        At ASL the archives are located in:
            /tr1/telemetry_days
          or
           /msd
        depending on how long since the event passed.
    """
    stime = eventtime - 3000.
    etime = eventtime + 3000. + lents

    # Current default is LH add BH later
    chan = 'LH'

    if (net in set(['GT'])) or (sta == 'KBL'):
        chan = 'BH'

    # Grab the data locations
    datalocs = ['/tr1/telemetry_days/', '/msd/']

    files = []
    # Grab the files and deal with edge cases
    for dataloc in datalocs:
        for year in range(stime.year, etime.year+1):
            for day in range(stime.julday, etime.julday+1):
                string = dataloc + net + '_' + sta + '/' + \
                        str(year) + '/' + '*' + \
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
                print('Unable to get data ' + curfile)
    st.merge(fill_value='latest')
    if debug:
        print('We have data')
    return st



def getcolor(chan, loc):
    """ Set the color of the trace in the plot depending on the channel
    or if it synthetic or observed data. """
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
    """ 
    calculate the correlation coefficient and lag time for the synthetic
    when compared to the observed data and write to a file.
    """
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
            print('No residual for' + tr.stats.station + ' ' + 'LH' + comp) 
    return


def getargs():
    """ Grab command line arguments to run synthetics. """
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

    
def procStream(st, inv, eventtime, tshift, freqmin, freqmax, corners, lents, hdur, debug=False):

    """ Deconvolve and filter each trace in the stream. """
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
# Here is where I (Adam) am mucking around
# Here is where Kim wonders if Adam is done mucking around.
            if synthetic:
                tr.data /= (10**9)
                tr.stats.starttime += float(tshift)
                win = signal.hann(int(2*hdur))
                tr.data = signal.convolve(tr.data,win, mode='same')/sum(win)
      
        except:
            print('Problem with response')
            st.remove(tr)
        try:    
            tr.taper(max_percentage=0.05, type='cosine')
            if not synthetic:
                tr.remove_response(inv)
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
    
def pltStream(stream, pltHandle, component,inv, cmtlat=None, cmtlon=None,
              minfre=None, maxfre=None, resDir=None, debug = False):
    """ Plot the synthetic and observed traces. """
    if component == 'Z':
        plt.subplot(3,1,1)
        title = stream[0].stats.network + ' ' + \
                stream[0].stats.station + ' '
        starttime = stream[0].stats.starttime
        stime = str(starttime.year) + ' ' + str(starttime.julday) + \
                ' ' + str(starttime.hour) + ':' + \
                str(starttime.minute) + ':' + \
                str("{0:.2f}".format(starttime.second))
        title += stime + ' '
        lat = inv[0][0].latitude
        lon = inv[0][0].longitude
        dist = gps2dist_azimuth(float(cmtlat), float(cmtlon), lat, lon)
        bazi = "{0:.1f}".format(dist[2])
        dist = "{0:.1f}".format(0.0089932*dist[0]/1000.)
        title += 'Dist:' + str(dist)
        title += ' BAzi:' + str(bazi) + ' ' 
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
# numpy.arange creates an evenly spaced array from 0, npts in the stream.
# dividing by sampling rate gives us a time.
    t = np.arange(0, stream[0].stats.npts)/stream[0].stats.sampling_rate
    stream.sort(['location'])
    for tr in stream.select(component=component):
        color = getcolor(tr.stats.channel, tr.stats.location)
        plt.plot(t, tr.data*(10**3), color, label = tr.stats.location +
                 ' ' + tr.stats.channel, linewidth=1)
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
    """ Return the sncl """
    nslc = (tr.id).split('.')
    return nslc[1], nslc[0], nslc[3], nslc[2]
                


# Start of the main part of the program
if __name__ == "__main__":
    
    # Eventually allow this to change in the arguments
    client = Client("IRIS")

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

    
    # Run through each of the event CMTs
    for synfile in parserval.syn:
     
        if synfile[-1] == '/':
            synfile = synfile[:-1]

        # Read in the CMT solution from the synthetic directory
        if debug:
            print("We are using local synthetics")
        if not os.path.isfile(synfile + '/CMTSOLUTION'):
            print("No CMT found")
            exit(0)
        try:
            cat = read_events(synfile + '/CMTSOLUTION')
        except:
            try:
                # Here we have a use case where we are missing a space
                copyfile(synfile + '/CMTSOLUTION', 'CMTTEMP')
                f=open('CMTTEMP','r')
                CMT = f.read()
                f.close()
                f = open('CMTTEMP', 'w')
                f.write(' ' + CMT)
                f.close()
                cat = read_events('CMTTEMP')
                os.remove('CMTTEMP')
            except:
                print("No CMT found")
                exit(0)
        if debug:
            print(cat)        
        
        cmtlat, cmtlon, eventtime, tshift, hdur = readcmt(cat)
        if debug:
            print(cat[0].origins)
        # Lets make a local results directory
        resultdir = parserval.resDir
        if resultdir[-1] == '/':
            resultdir = resultdir[:-1]

        if not os.path.exists(os.getcwd() + '/' + resultdir):
            os.mkdir(os.getcwd() + '/' + resultdir)
        evename = str(cat[0].resource_id).split('/')[-2]


        # Open a file to write the correlation statistics 
        statfile = open(os.getcwd() + '/' + resultdir + '/Results' + evename + net + '.csv' ,'w')
        statfile.write('net, sta, loc, chan, scalefac, lag, corr, time\n')

        if parserval.sta:
            if debug: 
                print("We are using a manual station list")
            stations = parserval.sta.split(",")
        else:
            stations = client.get_stations(network=net, starttime=eventtime, endtime=eventtime)
            stations = [ sta.code for sta in stations[0]]
        for sta in stations:
            #try:
            if True:
             
                st = getdata(net, sta, eventtime, lents)
                if len(st) == 0:
                    pass
                # Decimate any high sample rate traces to 1
                for tr in st:
                    if tr.stats.sampling_rate > 20:
                            tr.detrend('linear')
                            tr.taper(max_percentage=0.05, type='cosine')
                            tr.decimate(int(tr.stats.sampling_rate/4.))
                            tr.decimate(4)
                
                # We should get some metadata
                try:
                    inv = client.get_stations(starttime=eventtime, endtime=eventtime+5., network=net,
                                        sta = sta, channel="LH*", level="response")
                except:
                    continue
                # Lets go through each trace in the stream and deconvolve and filter
                st = procStream(st, inv, eventtime,0., userminfre, usermaxfre, filtercornerpoles, lents, hdur)
                # Lets check for reverse polarity and fix
                
                try:
                    st.rotate('->ZNE', inventory = inv)
                    st = choptocommon(st)
                except:
                    continue
                       
                for tr in st.select(channel="BH*"):
                    tr.stats.channel = tr.stats.channel.replace('B','L')

                files = glob.glob(synfile + '/' + sta + '.*.LX*.modes.sac')
                synstream = Stream()
                for curfile in files:
                   synstream += read(curfile)
                for tr in synstream:
                    tr.stats.channel = (tr.stats.channel).replace('LH','LX')
                    tr.stats.network = net  
                synstream = procStream(synstream, inv, eventtime, tshift, userminfre, usermaxfre, filtercornerpoles, lents, hdur)
                
                st += synstream

                
                # synplot is a plot handle that opens up a figure.
                synplot =  plt.figure(1)
                for comp in ["Z", "N", "E"]:
                    try:

                        pltStream(st, synplot, comp, inv,
                                  cmtlat=cmtlat, cmtlon=cmtlon,
                                  minfre=userminfre, maxfre=usermaxfre,
                                  resDir=resultdir)
                    except:
                        print('Problem with: ' + sta + ' plotting')
                    # Time to write some info into the statfile
                    try:
                        writestats(statfile, st, comp)
                    except:
                        print('Problem with: ' + sta)
        statfile.close()
