#!/usr/bin/perl -w

use strict ;

my $usage=<<USAGE;
Usage:

perl $0 <kog_fasta_blast_oufmt6_tab> <kog/whog> <COG/KOG>    >out.KOG.txt

Parameters :
	1-kog_fasta_blast_oufmt6_tab : blast to COG or KOG outfmt 6 output
	2-kog/whog : if cog , whog; if kog, kog;
	3-COG/KOG : which functional anno process you want to exec;


Download :
#fun.txt: `wget ftp://ftp.ncbi.nih.gov/pub/COG/KOG/fun.txt`
kog: `wget ftp://ftp.ncbi.nih.gov/pub/COG/KOG/kog`

ftp download advised by Thunder ;

Author: Changfu, Jia
Email : 1020160171\@qq.com
USAGE

die "$usage\n" unless @ARGV==3;

my $tab=shift;
my $kog=shift;
my $class=shift;

my (%cate, %KOG, %sp);


if ($class eq "KOG"){
$/="\n\n";
open KO, "$kog" or die "";
while (<KO>){
	chomp;
	my @line=split/\n/, $_;
	my $header= shift @line;
	$header =~ s/^\[//;
	my @header = split/] /, $header;
	#my $cate= $header[0];
	my $cate = substr($header[0],0,1);
	foreach my $l (@line){
		my @a=split/:\s+/, $l;
		$cate{$a[1]}=$cate;
		$KOG{$a[1]}=$header[1];
		$a[0]=~ s/^\s+//g;
		$sp{$a[1]}=$a[0];
		#print "$a[1]\t$a[0]\t$l\n";
	}
	#print "$cate\n";
}
close KO;
}elsif ($class eq "COG"){
$/="_______\n\n";
open CO, "$kog" or die "";
while (<CO>){
	chomp;
	my @line=split/\n/, $_;
	#print "${_}111\n";
	my $header= shift @line;
	#my $line2=join("\n", @line);
        $header =~ s/^\[//;
        my @header = split/] /, $header;

        my $cate = substr($header[0],0,1);
	#my @line2= split/:/, $line2;
	foreach my $l (@line){
		next unless $l =~ m/:/;
		my @a=split/:/, $l;
		#print "$a[1]\n";
		my @b=split/\s+/, $a[1];
		shift @b;
		foreach my $b (@b){
			#print "$a[0]\t$b\n";
               		$cate{$b}=$cate;
        	        $KOG{$b}=$header[1];
	                $sp{$b}=$a[0];

		}

	}


}
close CO;
}else {
	die "Third parameter must be \"COG\" or \"KOG\" \n $usage\n";
}



$/="\n";
open TAB, "$tab" or die "";
#open OUT, ">$tab.best" or die "";
my (%score, %gene);
while (<TAB>){
	chomp;
	my @line =split/\t/, $_;
	next unless $line[2] >40;
	next unless exists $sp{$line[1]};
	$line[0] =~ s/\|.*//g;
	$score{$line[0]}{$line[-1]}= $line[1];
}
close TAB;
#close OUT;

foreach my $gene (sort keys %score){
	my $c=1;
	foreach my $score (sort {$b <=> $a} keys %{$score{$gene}}){
		last unless $c==1 ;
		my $v=$score{$gene}{$score};
		$gene{$gene}=$v;
		#print "$gene\t$score\t$v\n";
		$c++;	
	}
}

#my (%cate, %KOG, %sp);

foreach my $gene (sort keys %gene){
	my $o=$gene{$gene};
	next unless exists $sp{$o};
	print "$gene\t$o\t$sp{$o}\t$cate{$o}\t$KOG{$o}\n";	
}







