#!/usr/bin/env python

######################
#UpdatesynCMT.py
#Adam Ringler
#
#This code takes the CMTmineos solutions and produces synthetics.  This
#code assumes you have an installation of MINEOS and it is setup with
#path information
#####################

import os
import glob
import sys
import shutil

debug = False

#Here is where your CMTs are as created by updateCMT.py
cmtdirpath = '/home/aringler/syncomppython'

#Here is where the MINEOS code resides
codepath = '/home/aringler/syncomppython'

def rungreen(event,codepath,currdir):
#This function runs the greens functions using the codepath and the event info

#Here is the length of the time series we are using
	tslen = str(8000)

#Lets setup the parameter file to run the Greens functions		
	parafile = open('parameter_file','w')
	parafile.write(codepath + '/auxfiles/longNEW\n')
	parafile.write(codepath + '/auxfiles/modefiles/db_list\n')
	parafile.write(event + '/currCMTmineos\n')

#Lets do the first 260 modes
	parafile.write('0 260\n')

	parafile.write(tslen + '\n')
	parafile.write('green\n')
	parafile.close()

#We can now feed the code into the greens program
	os.system(codepath + '/bin/green < parameter_file')
	os.remove('parameter_file')
	shutil.copy2(codepath + '/auxfiles/longNEW.site', currdir + '/green.site')			
	shutil.copy2(codepath + '/auxfiles/longNEW.sitechan', currdir + '/green.sitechan')
	return;


def runsyndat(event,codepath,currdir):
#This function runs the synthetics using syndat from Mineos
	if os.path.exists(currdir + 'Syndat.wfdisc'):
#Need to clean up the Syndat
			print 'Need to remove Syndat\n'
	parafile = open('parameter_file','w')
	parafile.write(event + '/currCMTmineos\n')
	parafile.write('0\n')
	parafile.write('green\n')
	parafile.write('Syndat\n')
	parafile.write('0\n')
	parafile.close()
	os.system(codepath + '/bin/syndat < parameter_file')
	os.remove(currdir + '/parameter_file')
	shutil.copy2(codepath + '/auxfiles/longNEW.site', currdir + '/Syndat.site')			
	shutil.copy2(codepath + '/auxfiles/longNEW.sitechan', currdir + '/Syndat.sitechan')
	os.system(codepath + '/bin/creat_origin ' + event + '/currCMTmineos Syndat')	
	return;







#Lets get all of the events for years starting with 1 and 2
eventlist = glob.glob(cmtdirpath + '/2*/*') + glob.glob(cmtdirpath + '/1*/*')
currdir = os.getcwd()

if debug:
	print eventlist
	print 'The current directory is:' + currdir

#Time to loop through all of the events and create synthetics for them
for ind, event in enumerate(eventlist):
	print 'On event ' + str(ind + 1) + ' of ' + str(len(eventlist)) 
	
	if debug:
		print event
	
#Okay we have an event if there are more than 2 files in the directory
#This should eventually be made more robust
	if len(os.listdir(event)) == 2:
		if debug:
			print 'No synthetics for this event'

#Run the greens functions		
		rungreen(event,codepath,currdir)

#Now we run the syndat piece of Mineos		
		runsyndat(event,codepath,currdir)
		
#Lets convert the synthetics to SAC format
		os.system(codepath + '/bin/cucss2sac Syndat Syns')
#Time to clean up stuff
		os.system('rm -r ' + currdir + '/Syndat.*')
		os.system('rm -r ' + currdir + '/green.*')
		synall = glob.glob(currdir + '/Syns/*.SAC')

#Here we are changing from H to X to get no network		
		for syncur in synall:
			os.rename(syncur, syncur.replace(' ',''))
			syncur=syncur.replace(' ','')
			syncurchan = syncur.replace(':','.')
			syncurchan = syncurchan.split('.')
			syncurchan[6] = syncurchan[6].replace('H','X')
			syncurchan = syncurchan[5] + '.XX.' + syncurchan[6] + '.modes.sac'
			if debug:
				print 'Old synthetic:' + syncur
				print 'New synthetic:' + syncurchan
				print 'Move to location:' + event + '/' + syncurchan
			os.system('mv ' + syncur + ' ' + event + '/' + syncurchan)
		os.system('rm -r ' + currdir + '/Syns')

#Here we setup our process_syn.pl code to do the pre-processing
#This should be removed and done in python later as it is really sloppy
		synprostr = '-S -m ' + event + '/CMTSOLUTION '
		synprostr = synprostr + '-s 1.0 -l 0/4000 -t 40/400 -x proc ' 
		synprostr = synprostr + event + '/*modes.sac '	
		print synprostr
		os.system(codepath + '/bin/process_syn.pl ' + synprostr )
		os.system('rm -r ' + event + '/*modes.sac')
		os.system('chmod -R 755 ' + event)
sys.exit(0)
	



