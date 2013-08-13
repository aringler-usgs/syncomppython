#!/usr/bin/perl

#This script converts the date to the jday


if($#ARGV != 2){
	print "Usage Year Month Day\n";
	exit;
}

$year = $ARGV[0];
$mon = $ARGV[1];
$day = $ARGV[2];



if($year % 400 == 0){
	$lyear = "TRUE";
}
elsif($year % 100 == 0){
       $lyear = "FALSE";
}
elsif($year % 4 == 0){
       $lyear = "TRUE";
}
else{
       $lyear = "FALSE";
}


if($lyear eq "TRUE"){
	if($mon == 1){
		$jday = $day;
	}
	elsif($mon == 2){
		$jday = $day + 31;
	}
	elsif($mon == 3){	
		$jday = $day + 60;
	}
	elsif($mon == 4){
		$jday = $day + 91;
	}
	elsif($mon == 5){	
		$jday = $day + 121;
	}
	elsif($mon == 6){
		$jday = $day + 152;
	}
	elsif($mon == 7){	
		$jday = $day + 182;
	}
	elsif($mon == 8){
		$jday = $day + 213;
	}
	elsif($mon == 9){	
		$jday = $day + 244;
	}
	elsif($mon == 10){	
		$jday = $day + 274;
	}
	elsif($mon == 11){
		$jday = $day + 305;
	}
	elsif($mon == 12){	
		$jday = $day + 335;
	}

}
else{

	if($mon == 1){
		$jday = $day;
	}
	elsif($mon == 2){
		$jday = $day + 31;
	}
	elsif($mon == 3){	
		$jday = $day + 59;
	}
	elsif($mon == 4){
		$jday = $day + 90;
	}
	elsif($mon == 5){	
		$jday = $day + 120;
	}
	elsif($mon == 6){
		$jday = $day + 151;
	}
	elsif($mon == 7){	
		$jday = $day + 181;
	}
	elsif($mon == 8){
		$jday = $day + 212;
	}
	elsif($mon == 9){	
		$jday = $day + 243;
	}
	elsif($mon == 10){	
		$jday = $day + 273;
	}
	elsif($mon == 11){
		$jday = $day + 304;
	}
	elsif($mon == 12){	
		$jday = $day + 334;
	}










}





print "${jday}\n";
