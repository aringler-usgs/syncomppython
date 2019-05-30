#! /usr/bin/env python

"""
Compute synthetic seismograms.

This script has many directory dependencies that may require updating of
your .bashrc $PATH variable.  You need to have /usr/local/bin in order to 
run the mineos functions.  Try a "which green" at a command prompt to make sure that it is in your user path and/or exists on the machine you are running.

There are also precomputed mode files which you will need either a copy or a 
symbolic link in your codepath. Make sure that you update the codepath variable 
to reflect your directory structure. At the moment it will look for the infor-
mation in $HOME/synInfo.  synInfo should have subdirectories:
     auxfiles/ - contains station information files and modefiles/ directory.
    modefiles/ - contains the normal mode information for calculating the
                 synthetic seismogram.
"""
import os
import glob
import sys
import shutil

debug = True

# this is where the program expects to find the mineos input file 
# containing the CMT event information
cmtdirpath = '/home/aringler/data_stuff/syncomppython/SYNTHETICS'
# this is where the program expects to find the station information and the
# mode information.
codepath = '/home/aringler/data_stuff/syncomppython/synInfo'
eventlist = glob.glob(cmtdirpath + '/2*/*') + glob.glob(cmtdirpath + '/1*/*')
currdir = os.getcwd()
ind=0;
if debug:
    print eventlist
    print 'The current directory is:' + currdir
for event in eventlist:
    ind += 1
    print 'On event ' + str(ind) + ' of ' + str(len(eventlist)) 
    if debug:
        print event
    if len(os.listdir(event)) == 2:
        if debug:
            print 'No synthetics for this event'
# Set up the input file to run the mineos greens function.        
# the output from the greens function is a *.wf_disc file used in 
# syndat
# inputs needed are precomputed eigenfunctions from mineos_bran
# and eigcon.  
        parafile = open('parameter_file','w')
        if debug:
            print 'writing parameters'
# this is the path to the database files (directory containing .site and
# .sitechan files)
        parafile.write(codepath + '/auxfiles/longNEW\n')
# this is the file within the normal modes database
        parafile.write(codepath + '/auxfiles/modefiles/db_list\n')
# this is the file with the CMT information
        parafile.write(event + '/CMTSOLUTIONmineos\n')
# this is the minumum and maximum frequencies
        parafile.write('0. 260.\n')
# this is the number of pts in the greens function.  This must be >= 30K
        parafile.write('8000\n')
# this is the green functions output file name
        parafile.write('green\n')
        parafile.close()
# this will take the parameter file created above and calculate the 
# green's functions
        os.system('bin/green < parameter_file')
# then we remove the file.

        os.remove('parameter_file')

# this is the path to the station information.
        shutil.copy2(codepath + '/auxfiles/longNEW.site', currdir + '/green.site')            
        shutil.copy2(codepath + '/auxfiles/longNEW.sitechan', currdir + '/green.sitechan')
        
# Check to see if we need to clean up the Syndat
        if os.path.exists(currdir + 'Syndat.wfdisc'):
            print 'Need to remove Syndat\n'
# build the input parameter file for cucss2sac - transforms file to sac format
        parafile = open('parameter_file','w')
        parafile.write(event + '/CMTSOLUTIONmineos\n')
        parafile.write('0\n')
        parafile.write('green\n')
        parafile.write('Syndat\n')
        parafile.write('0\n')
        parafile.close()
        os.system('bin/syndat < parameter_file')
        os.remove(currdir + '/parameter_file')
        shutil.copy2(codepath + '/auxfiles/longNEW.site', currdir + '/Syndat.site')            
        shutil.copy2(codepath + '/auxfiles/longNEW.sitechan', currdir + '/Syndat.sitechan')
        os.system('bin/creat_origin ' + event + '/CMTSOLUTIONmineos Syndat')

        os.system('bin/cucss2sac Syndat Syns')
#Time to clean up stuff
        os.system('rm -r ' + currdir + '/Syndat.*')
        os.system('rm -r ' + currdir + '/green.*')
#rename the output files to something we like better
#first create a list
        synall = glob.glob(currdir + '/Syns/*.SAC')
        for syncur in synall:
# get rid of spaces and replace with 0 - makes all of doy's 3 chars
            os.rename(syncur, syncur.replace(' ','0'))
            syncur=syncur.replace(' ','0')
# replace colons with periods
            syncurchan = syncur.replace(':','.')
            syncurchan = syncurchan.split('.')
# rename the channels, using X for a synthetic
            syncurchan[6] = syncurchan[6].replace('H','X')
# add a different ending
            syncurchan = syncurchan[5] + '.XX.' + syncurchan[6] + '.modes.sac'
            if debug:
                print 'Old synthetic:' + syncur
                print 'New synthetic:' + syncurchan
                print 'Move to location:' + event + '/' + syncurchan
#  move things around 
            os.system('cp ' + syncur + ' ' + event + '/' + syncurchan)
# do some clean up
        sys.exit()        
        os.system('rm -r ' + currdir + '/Syns')
# add some info into the synthetic file headers using a perl script. if
# you are curious about these options take a look at the script.
        synprostr = '-S -m ' + event + '/CMTSOLUTION '
        synprostr = synprostr + '-s 1.0 -l 0/4000 -t 40/400 -x proc ' 
        synprostr = synprostr + event + '/*modes.sac '    
        print synprostr
## make sure the process_syn.pl is in your path
        os.system('bin/process_syn.pl ' + synprostr )
        os.system('rm -r ' + event + '/*modes.sac')
        os.system('chmod -R 755 ' + event)
        sys.exit()
