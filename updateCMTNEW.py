#!usr/bin/env python
import os
import sys
import datetime
import math
from obspy import read_events
import obspy
from obspy.io.xseed import Parser
from obspy.core import UTCDateTime
import shutil
import glob

debug = True

def cmt_to_mineos(filename, debug = True):
    try:
        cat = obspy.read_events(filename)
    except:
        with open(filename, 'r+') as f:
            old = f.read()
            f.seek(0)
            f.write(" " + old)
        cat = obspy.read_events(filename)
    if debug:
        print(cat[0].origins)
    fcmt = open(filename)
    line= fcmt.readlines()
    tshift = line[2].replace('time shift:','')
    tshift = float(tshift.strip())
    cat[0].origins[0].time += -tshift
    f = open(filename + 'mineos', 'w')
    eveid = str(cat[0].origins[0].time.year)[2:] 
    eveid += str(cat[0].origins[0].time.month).zfill(2)
    eveid += str(cat[0].origins[0].time.day).zfill(2)
    eveid += str(cat[0].origins[0].time.hour).zfill(2)
    #eveid += str(cat[0].origins[0].time.minute).zfill(2)
    #eveid += 'A'
    f.write(eveid + ' ')
    time = str(cat[0].origins[0].time.year) + ' ' 
    time += str(cat[0].origins[0].time.julday) + ' '
    time += str(cat[0].origins[0].time.hour).zfill(2) + ' '
    time += str(cat[0].origins[0].time.minute).zfill(2) + ' '
    time += str(cat[0].origins[0].time.second).zfill(2) + '.00 '
    f.write(time)
    lat = str('{:.2f}'.format(cat[0].origins[0].latitude)).ljust(5)
    f.write(lat)
    lon = str(cat[0].origins[0].longitude).rjust(8)
    f.write(lon)
    depth = ' ' + str('{:.1f}'.format(cat[0].origins[1].depth/1000.))
    if debug:
        print('Here is the depth:' + depth)
    f.write(depth)
    # Now write the time step in seconds
    f.write(' 20.0')
    if debug:
        print(cat[0].focal_mechanisms[0])
    mom = cat[0].focal_mechanisms[0].moment_tensor
    f.write(' ' + str('{:.1f}'.format(mom.source_time_function['duration']/2.)))
    f.write(' ' + str('{:.3e}'.format(mom.scalar_moment*10**7)).replace('+',''))
    exp = int(str(mom.scalar_moment).split('+')[1])
    f.write(' ' + str(mom.tensor['m_rr']/(10**exp)))
    f.write(' ' + str(mom.tensor['m_tt']/(10**exp)))
    f.write(' ' + str(mom.tensor['m_pp']/(10**exp)))
    f.write(' ' + str(mom.tensor['m_rt']/(10**exp)))
    f.write(' ' + str(mom.tensor['m_rp']/(10**exp)))
    f.write(' ' + str(mom.tensor['m_tp']/(10**exp)))
    f.write(' ' + str(1.0) + 'e' + str(exp + 7))
    a = obspy.imaging.beachball.MomentTensor(mom.tensor['m_rr'], mom.tensor['m_tt'], 
    mom.tensor['m_pp'], mom.tensor['m_rt'], mom.tensor['m_rp'], mom.tensor['m_tp'], 1.)
    
    fp = obspy.imaging.beachball.mt2plane(a)
    f.write(' ' + str(int(round(fp.strike))))
    f.write(' ' + str(int(round(fp.dip))))
    f.write(' ' + str(int(round(fp.rake-360.))))
    if debug:
        print(fp.strike)
    ax = obspy.imaging.beachball.aux_plane(fp.strike, fp.dip, fp.rake)
    f.write(' ' + str(int(round(ax[0]))))
    f.write(' ' + str(int(round(ax[1]))))
    f.write(' ' + str(int(round(ax[2]))) + '\n')
        
    f.close()
    return
    




# this is the path where the synthetics will be created.
cmtdirpath = '/home/aringler/src/fix_syncomp/syncomppython/SYNTHETICS'
# this is the path where it will output information about all of the data it is 
# getting. at the moment it is set to a directory in the user's home
codepath = '/home/aringler/src/fix_syncomp/syncomppython/synInfo'
minmag= 7.0

#Download the latest CMT files
currdir = os.getcwd()
if not os.path.exists(codepath + '/croncode'):
    os.mkdir(codepath + '/croncode')
os.chdir(codepath + '/croncode')
os.system('wget -N http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk')
os.chdir(currdir)

#Read in the CMT file

cat = read_events(codepath + '/croncode/qcmt.ndk')
for eve in cat:
    if debug:
        print(eve)
    eveid = eve.resource_id.id.split('/')[2]
    evemag = eve.magnitudes[0].mag
    etime = eve.origins[0].time
    print(etime)
    # Flow control break to next if too small
    if evemag < minmag:
        continue
    
    
    if not os.path.exists(cmtdirpath + '/' + str(etime.year)):
        os.mkdir(cmtdirpath + '/' + str(etime.year))
    cmtdire = cmtdirpath + '/' + str(etime.year) + '/' + eveid
    cmtdire = cmtdire.strip()
    if not os.path.exists(cmtdire):
        os.mkdir(cmtdire)
    eve.write(cmtdire + '/CMTSOLUTION', format='CMTSOLUTION')
    cmt_to_mineos(cmtdire + '/CMTSOLUTION')
