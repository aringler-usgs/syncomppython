Synthetic Comparison Python
=============

Synthetics comparison code

This code compares MINEOS synthetics to data and plots the results using obspy.  The code is mostly written for the internal data directory structure of the Albuquerque Seismological Laboratory.


Usage
=============

The code uses and argument line parser to change various parameters in the comparison.

./syncomp.py -h
usage: syncomp.py [-h] -n NETWORK -resDir RESDIR -syn SYN [SYN ...] [-sta STA]
                  [-tslen LENTS] [-debug] [-dataloc]
                  [-filter FILTER FILTER FILTER]

Program to compare long-period event synthetics to data

optional arguments:

  -h, --help            show this help message and exit

  -n NETWORK            Network name Example: IU

  -resDir RESDIR        Result directory name Example: blah

  -syn SYN [SYN ...]    Synthetics directory location Example:
                        /SYNTHETICS/2014/C201401*

  -sta STA              Stations to use Example with a comma (,) separator :
                        TUC,ANMO

  -tslen LENTS          Length of time series in seconds Example: 4000,
                        default is 4000 s

  -debug                Run in debug mode

  -dataloc              Use /xs0 data location, otherwise use /tr1 also

  -filter FILTER FILTER FILTER
                        Filter parameters using minimum period maximum period
                        and number of corners Example: 100 200 4, default is
                        100 400 2



Example 
=============

./syncomp.py -n CU -resDir blah3 -syn /SYNTHETICS/2014/C2014030* -sta GRTK -filter 100 300 4
