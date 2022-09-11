#!/usr/bin/perl -w

use strict;


my $domain=shift;
my $qry=shift;
my $ref=shift;
my $hal=shift;



my $mark=1;

open COM, ">${qry}_${ref}/command/${qry}_qry_$ref.liftover" or die "";

open IN, "$domain" or die "";
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	my $boudA1=$line[1]-10000;
	my $boudA2=$line[1]+10000;
	my $boudB1=$line[2]-10000;
	my $boudB2=$line[2]+10000;
	open DOM, ">${qry}_${ref}/$qry.TAD_domain.$mark.bed";
	print DOM "$_\n";
	close DOM;
	open BOA, ">${qry}_${ref}/$qry.TAD_boundaryA.$mark.bed";
	print BOA "$line[0]\t$boudA1\t$boudA2\n";
	close BOA;
	open BOB, ">${qry}_${ref}/$qry.TAD_boundaryB.$mark.bed";
	print BOB "$line[0]\t$boudB1\t$boudB2\n";
	close BOB;
	print COM "halLiftover $hal $qry ${qry}_${ref}/$qry.TAD_domain.$mark.bed  $ref ${qry}_${ref}/$qry.TAD_domain.$mark.liftover\nhalLiftover $hal $qry ${qry}_${ref}/$qry.TAD_boundaryA.$mark.bed $ref ${qry}_${ref}/$qry.TAD_boundaryA.$mark.liftover \n halLiftover $hal $qry ${qry}_${ref}/$qry.TAD_boundaryB.$mark.bed $ref ${qry}_${ref}/$qry.TAD_boundaryB.$mark.liftover \n" ;
	$mark++;
}




