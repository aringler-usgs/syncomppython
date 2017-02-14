#!/usr/bin/env python
import glob
import os

debug = True


#Here is where you can have the results placed.  You will need to change this
resdir = os.getenv('HOME')+'/synthetics/results'

#Here is where the synthetics reside.  You probably do not need to change this
syndir = '/SYNTHETICS'

#Here is the path to where syncomp.py resides this will be in your local git repostiry
# you may need to update this depending on where you keep your scripts
syncompdir = os.getenv('HOME')+'/python/syncomppython'
nets = ['IU','IC','CU','US','IW','NE','II']

if debug:
	syndirs = glob.glob('/SYNTHETICS/*/*')
else:
	syndirs = glob.glob('SYNTHETICS/*/*')

dondirs = glob.glob(resdir)

if debug:
	print syndirs
	print dondirs

curdir = os.getcwd()


for event in syndirs:
	dirinfo = event.split('/')
	if not os.path.exists(resdir + '/' + dirinfo[2]):
		if debug:
			print 'need to make a directory'
		os.system('mkdir ' + resdir + '/' + dirinfo[2])
	chgdir = resdir + '/' + dirinfo[2]
	os.chdir(chgdir)
	if debug:
		print 'New directory: ' + chgdir
	for curnet in nets:
		if debug:
			print 'Current network: ' + curnet
		chkpath = resdir + '/' + dirinfo[2] + '/' + \
			dirinfo[3] + '/Results' + dirinfo[3] + curnet + '.csv'
		if debug:
			print 'Here is the path we are checking: ' + chkpath
		if not os.path.exists(chkpath):
			runsyn = syncompdir + '/syncomp.py ' + syndir + '/' + dirinfo[2] + \
				'/' + dirinfo[3] + ' ' + dirinfo[3] + ' ' + curnet
			if debug:
				print 'Running: ' + runsyn
			os.system(runsyn)






os.chdir(curdir)
