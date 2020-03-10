#!/usr/bin/perl

#Removes 4nt randomers from each side of sequencing reads in fastq files
#Appends randomers to read name
#Takes one argument, the file name
#Written by Emanuel Wyler, emanuel.wyler@alumni.ethz.ch

use strict;

#print "Hello World\n";

my $filename =$ARGV[0];
my $N1;
my $N2;
my $filenameout = $filename;
$filenameout =~ s/(.*)\.([^\.]*)$/$1-2xNNNN\.fastq/;
#my $filenameout =$ARGV[1];

#print ("\nFilename in:",$filename,"\n");

#print ("\nFilename out: ",$filenameout,"\n");

open(FILE,$filename); 
open(FOUT,">$filenameout");

while (my $line1 = <FILE>) {
	my $line2 = <FILE>;
	my $line3 = <FILE>;
	my $line4 = <FILE>;
	chomp $line1;
	chomp $line2;
	chomp $line3;
	chomp $line4;
	$N1=substr($line2, 0, 4);
	$N2=substr($line2, length($line2)-4, 4);
	$line2 = substr($line2, 4, length($line2)-8);
	$line4 = substr($line4, 4, length($line4)-8);
	$line1 = $line1 . "_" . "$N1:$N2";
	print FOUT "$line1\n$line2\n$line3\n$line4\n";
}

close FILE;
close FOUT;
