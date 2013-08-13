#!/usr/bin/env python

import os
import glob
import sys
import shutil

debug = bool(1)

cmtdirpath = '/home/aringler/syncomppython'
codepath = '/home/aringler/syncomppython'
eventlist = glob.glob(cmtdirpath + '/2*/*') + glob.glob(cmtdirpath + '/1*/*')
currdir = os.getcwd()
ind=0;
if debug == bool(1):
	print eventlist
	print 'The current directory is:' + currdir
for event in eventlist:
	ind += 1
	print 'On event ' + str(ind) + ' of ' + str(len(eventlist)) 
	if debug == bool(1):
		print event
	if len(os.listdir(event)) == 2:
		if debug == bool(1):
			print 'No synthetics for this event'
		
		parafile = open('parameter_file','w')
		parafile.write(codepath + '/auxfiles/longNEW\n')
		parafile.write(codepath + '/auxfiles/modefiles/db_list\n')
		parafile.write(event + '/currCMTmineos\n')
		parafile.write('0 260\n')
		parafile.write('8000\n')
		parafile.write('green\n')
		parafile.close()
		os.system(codepath + '/bin/green < parameter_file')
		os.remove('parameter_file')
		shutil.copy2(codepath + '/auxfiles/longNEW.site', currdir + '/green.site')			
		shutil.copy2(codepath + '/auxfiles/longNEW.sitechan', currdir + '/green.sitechan')
		
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
		
		os.system(codepath + '/bin/cucss2sac Syndat Syns')
#Time to clean up stuff
		os.system('rm -r ' + currdir + '/Syndat.*')
		os.system('rm -r ' + currdir + '/green.*')
		synall = glob.glob(currdir + '/Syns/*.SAC')
		for syncur in synall:
			os.rename(syncur, syncur.replace(' ',''))
			syncur=syncur.replace(' ','')
			syncurchan = syncur.replace(':','.')
			syncurchan = syncurchan.split('.')
			syncurchan[6] = syncurchan[6].replace('H','X')
			syncurchan = syncurchan[5] + '.XX.' + syncurchan[6] + '.modes.sac'
			if debug == bool(1):
				print 'Old synthetic:' + syncur
				print 'New synthetic:' + syncurchan
				print 'Move to location:' + event + '/' + syncurchan
			os.system('mv ' + syncur + ' ' + event + '/' + syncurchan)
		os.system('rm -r ' + currdir + '/Syns')
		synprostr = '-S -m ' + event + '/CMTSOLUTION '
		synprostr = synprostr + '-s 1.0 -l 0/4000 -t 40/400 -x proc ' 
		synprostr = synprostr + event + '/*modes.sac '	
		print synprostr
		os.system(codepath + '/bin/process_syn.pl ' + synprostr )
		os.system('rm -r ' + event + '/*modes.sac')
		os.system('chmod -R 755 ' + event)
sys.exit(0)
	



