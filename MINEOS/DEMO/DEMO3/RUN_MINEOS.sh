#!/bin/sh
#
# Example of the synthetic seismogram's generation script.
#
# Usage: RUN_MINEOS.sh model_name
# 
# Available model names for DEMO version:
# prem_noocean  prem_ocean CPacific NRussia
#
#=========================================================
#
# function creat_orig
#
creat_origin() {
time=`awk '{if(NR == 1){yd=$2-1970;vis=int((yd+1)/4); s=(365*yd+vis+$3-1)*86400.0+($4*60+$5)*60+$6; printf("%17.5f", s);}}' < $1`
awk 'BEGIN{t='"$time"';}{if(NR == 1) \
printf("%9.4f %9.4f %9.4f %17.5f %8d %8d %8d %4d %4d %4d %8d %8d %-7.7s %9.4f %-1.1s %7.2f %8d %7.2f %8d %7.2f %8d %-15.15s %-15.14s %8d %-17.17s\n", \
$7,$8,$9,t,1,-1,$2*1000+$3,-1,-1,-1,-1,-1,"-",-999.0000,"-",-999.0,-1,-999.0,-1, \
-999.00,-1,"-","PDE & Hvd CMT",-1,-1); \
}' $1 > $2.origin ;
}
#
#   Main procedure
#
if test "$#" != 1; then
echo " Usage: RUN_MINEOS.sh model_name"
exit
fi
model=$1                # setup 1-D model name
# check model name
flg=0
for f in \
CPacific NRussia prem_noocean prem_noocean_na prem_ocean prem_noocean_1ln.txt
do
if test "$f" = $1; then
flg=1
fi
done
if test "$flg" = 0; then
echo "Model name $1 is wrong, allowed names are:"
echo "CPacific NRussia prem_noocean prem_noocean_na prem_ocean prem_noocean_1ln.txt"
exit
fi
#=========================================================
# 1. run minos_bran program for fundamental S  mode,
# where,  n=0, 0 < f <  0.2 Hz,
#
echo "Step 1:  minos_bran runs for S modes ....................."
echo "============== Program minos_bran =================="
null=`ls * | grep '_S$'`
if test "X$null" != X; then 
rm -f $null
fi
time minos_bran << EOF
../models/${model}.txt
${model}_S
e${model}_S
1.0e-10 1
3
2 8000 0.0 200.0 0 0
EOF
#=========================================================
# 2. run minos_bran program for fundamental T  mode,
# where,  n=0, 0 < f <  0.2 Hz,
#
echo "Step 2: minos_bran runs for T modes ....................."
echo "============== Program minos_bran =================="
null=`ls * | grep '_T$'`
if test "X$null" != X; then 
rm -f $null
fi
time minos_bran << EOF
../models/${model}.txt
${model}_T
e${model}_T
1.0e-10 1
2
2 8000 0.0 200.0 0 0
EOF
#============================================================
# 3. Convert minos_bran results to .eigen relation (S mode)
#
echo "Step 3: eigen for S ....................................."
if test -f test_S.eigen; then
 rm -rf test_S.*
fi
time eigcon << EOF
3
../models/${model}.txt
1000
${model}_S
e${model}_S
test_S
EOF
#============================================================
# 4. Convert minos_bran results to .eigen relation (T mode)
echo "Step 4: eigen for T ....................................."
if test -f test_T.eigen; then
 rm -rf test_T.*
fi
time eigcon << EOF
2
../models/${model}.txt
1000
${model}_T
e${model}_T
test_T
EOF
#=========================================================
# 5. Evaluate green functions for given sitechan relation
echo "Step 5: green functions evaluation ........................."
if test -f green.wfdisc; then
 rm -rf green.*
fi
time green << EOF
short
db_list
china_cmt_event
10 260
8000
green
EOF
cp -p short.site green.site
cp -p short.sitechan green.sitechan
# create origin relation for data base green
creat_origin china_cmt_event green
#============================================================
# 6. Synthetic data construction
echo "Step 6: synthetic seismogram construction .................."
if test -f Syndat.wfdisc; then
 rm -rf Syndat.*
fi
time syndat << EOF
china_cmt_event
0
green
Syndat
0
EOF
cp -p short.site Syndat.site
cp -p short.sitechan Syndat.sitechan
creat_origin china_cmt_event Syndat
exit
