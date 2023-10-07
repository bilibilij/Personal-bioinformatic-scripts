#!/usr/bin/perl -w

while (<>){
	chomp;
	print "$1\n" if m/>(\S+)/;
}
