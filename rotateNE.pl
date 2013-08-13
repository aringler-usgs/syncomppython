#!/usr/bin/perl

if($#ARGV != 2){
	print "Usage: Z file 1 file 2 file";
	exit;
}


$compZ = $ARGV[0];
$comp1 = $ARGV[1];
$comp2 = $ARGV[2];


$sta = $compZ;
@sta = split(/\./,$sta);
$sta = @sta[0];


open(SAC,"|gsac ") || die("Can't open sac\n");

print SAC "r $compZ $comp1 $comp2\n";
print SAC "rot3 to 0\n";
print SAC "w\n";
print SAC "echo off\n quit\n";
close(SAC);

$comp1 =~ s/LH1/LHN/g;
$comp2 =~ s/LH2/LHE/g;

`mv ${sta}LH000 $comp1`;
`mv ${sta}LH090 $comp2`;
`mv ${sta}LHZ $compZ`;



