#!/usr/bin/perl
sub ltrim($);
sub rtrim($);

if($#ARGV  != 1) {
	print "Usage: Syntheticlocation ResultsName\n";
	exit;
}
$debug = "TRUE";
$synfile = $ARGV[0];
$resultdir = $ARGV[1];

if($debug eq "TRUE"){	
	$temp = substr($synfile,-6,6);
	print "${temp}\n";
}	

if(substr($synfile,-6,6) eq "tar.gz"){
	if($debug eq "TRUE"){
		print "Extract Princeton Synthetics\n";
	}
	`tar -xzf $synfile`;
	@synpiece = split(/\//,$synfile);
	$synpiece = @synpiece[length(@synpiece) - 1];
	$synpiece = substr($synpiece,0,-11);
	unless(-d $synpiece){	
		`mkdir $synpiece`;
	}
	`mv *.modes.sac $synpiece`;
	`mv CMTSOLUTION $synpiece`;
	$synfile = $synpiece;
	if(-f "STATIONS"){
		`rm STATIONS`;
	}
}



open stafile, "<stationlist"; 
while(<stafile>){
	$stainf = $_;
	$stainf =~ s/\n//g;
	$stainf=rtrim(ltrim($stainf));
	rtrim($stainf);
	if($debug eq "TRUE"){
		print "$stainf \n";
	}
	unless(-d $resultdir){
		`mkdir ${resultdir}`;
	}
	$net = substr($stainf,0,2);
	$sta = substr($stainf,3,length($stainf));
	if($debug eq "TRUE"){
		print "Station: $sta Network: $net\n";
	}
	
	unless(-f "${synfile}/CMTSOLUTION"){
		print "Can't find CMTSOLUTION\n";
		exit;
	}
	`cp ${synfile}/CMTSOLUTION .`;
	open(CMT,"CMTSOLUTION");
	$time = <CMT>;
	$year = substr($time,4,4); 
	$mon = substr($time,9,2);
	$day = substr($time,12,2);
	$mon = ltrim($mon);
	$day = ltrim($day);
	if($mon < 10){
		$mon = "0${mon}";
	} 

	if($day < 10){
		$day = "0${day}";
	}
	$jday = `date2jday.pl $year $mon $day`;
	chomp($jday);
	$jdayp1 = $jday + 1;
	if($jday < 100){
		$jday = "0${jday}";
	}
	if($jdayp1 < 100 && $jdayp1 > 10){
		$jdayp1 = "0${jdayp1}";
	}
	elsif($jdayp1 < 10){
		$jdayp1 = "00${jdayp1}";
	}
		
	
	@syns = glob("${synfile}/*${sta}*");
	if($debug eq "TRUE"){
		print 	"${synfile}/*${sta}*\n"
	}
	foreach(@syns){
		`cp $_ .`;
	}
	if($net eq "US" || $net eq "IW" || $net eq "NE"){
		`cat /xs1/seed/${net}_${sta}/${year}/${year}_${jday}*/*LH*.seed > ${sta}${net}.seed`;
		`cat /xs1/seed/${net}_${sta}/${year}/${year}_${jdayp1}*/*LH*.seed >> ${sta}${net}.seed`;
		`rdseed -f ${sta}${net}.seed -g /APPS/metadata/SEED/${net}.dataless -d -o 1`;
		`rm ${sta}${net}.seed`;
		
	}
	elsif($net eq "IU" || $net eq "CU" || $net eq "IC"){
		`cat /xs0/seed/${net}_${sta}/${year}/${year}_${jday}*/*LH*.seed > ${sta}${net}.seed`;
		`cat /xs0/seed/${net}_${sta}/${year}/${year}_${jdayp1}*/*LH*.seed >> ${sta}${net}.seed`;
		`rdseed -f ${sta}${net}.seed -g /APPS/metadata/SEED/_GSN.dataless -d -o 1`;
		`rm ${sta}${net}.seed`;
	}
	elsif($net eq "II"){
		`cat /TEST_ARCHIVE/${net}_${sta}/${year}/${year}_${jday}*/*LH*.seed > ${sta}${net}.seed`;
		`cat /TEST_ARCHIVE/${net}_${sta}/${year}/${year}_${jdayp1}*/*LH*.seed >> ${sta}${net}.seed`;
		`rdseed -f ${sta}${net}.seed -g /APPS/metadata/SEED/_GSN.dataless -d -o 1`;
		`rm ${sta}${net}.seed`;
	}
	@rdseederr = glob("rdseed.err*");
	foreach(@rdseederr){
		`rm $_`;
	}
	@datafiles = glob("*${sta}*Q.SAC *${sta}*D.SAC");
	foreach (@datafiles){
		$cursac = $_;
		chomp($cursac);
		@cursaccomps = split(/\./,$cursac);
		push(@chans, @cursaccomps[9]);
		if(length(@cursaccomps[8]) < 1){
			push(@locs, "NA");
		}	
		else{
			push(@locs, @cursaccomps[8]);	
		}
	}
	my %seen;
	my @uniquelocs = grep { ! $seen{$_}++ } @locs;
	my %seen;
	my @uniquechans = grep { ! $seen{$_}++ } @chans;
	

	foreach (@uniquelocs){
		$curloc = $_;
		if($debug eq "TRUE"){
			print "Here is the location: ${curloc}\n";
		}
		foreach(@uniquechans){	
			$curchan = $_;
			open(SAC,"|gsac") || die("Can't open sac.\n");
			if($curloc eq "NA"){
				print SAC "r *${sta}*${curchan}*.SAC\n";
			}	
			else{
				print SAC "r *${sta}*${curloc}*${curchan}*.SAC\n";
			}
			print SAC "merge\n";
			print SAC "w\n";
			print SAC "echo off\n quit\n";
			close(SAC);
			if($curloc eq "NA"){
				`mv ${sta}${curchan} ${sta}.${net}..${curchan}.SAC`;
				`/APPS/dcc/bin/getSACresploc ${net} ${sta} NA ${curchan} ${year} ${jday}`;
			}
			else{
				`mv ${sta}${curchan} ${sta}.${net}.${curloc}.${curchan}.SAC`;
				`/APPS/dcc/bin/getSACresploc ${net} ${sta} ${curloc} ${curchan} ${year} ${jday}`;
			}
		}
		
	}
	@olddata = glob("*${sta}*Q.SAC *${sta}*D.SAC");
	foreach(@olddata){
		`rm $_`;
	}
	`./process_data.pl -m CMTSOLUTION -l 1/3999 -i -p -t 40/400 -x proc *.SAC`;
	@syndata = glob("*.${sta}.*.modes.sac");
	`./process_syn.pl -m CMTSOLUTION -l 0/4000 -t 40/400 -S -x proc *.sac`;
	foreach(@syndata){		
		$cursyn = $_;		
		`rm $cursyn`;
		@tempchan = split(/\./,$cursyn);
		$tempchan = @tempchan[2];
		`mv ${tempchan}.proc ${sta}.XX.modes.sac.proc`;
	}



	@olddata = glob("*${sta}*.SAC");
	foreach(@olddata){
		`rm $_`;
	}
	@sacpzs = glob("SAC_PZs*${sta}*");
	foreach(@sacpzs){
		`rm $_`;
	}
	@hor12comps = glob("${sta}.${net}.*.LH1.SAC.proc");
	foreach(@hor12comps){
		$hor1 = $_;
		$hor2 = $hor1;
		$hor2 =~ s/LH1/LH2/g;
		$horZ = $hor1;
		$horZ =~ s/LH1/LHZ/g;
		`./rotateNE.pl $horZ $hor1 $hor2`;
		`rm $hor1 $hor2`;
	}
	`mv *${sta}*.SAC.proc $resultdir`;
	`mv *${sta}*.sac.proc $resultdir`;
}

if(-f "tmp.log"){
	`rm tmp.log`;
}
if(-f "ttimes.lst"){
	`rm ttimes.lst`;
}
`mv CMTSOLUTION $resultdir`;
if(-e $synpiece){
	`rm -r $synpiece`;
}

if(-f "currCMTmineos"){
	`rm currCMTmineos`;
}



# Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}
