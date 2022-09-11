#!/usr/bin/envs perl -w

use strict;

my $usage=<<USAGE;
perl $0 nlr.putative.fasta nlrparser3.out.txt > intact.fa

USAGE

die "$usage\n" if  @ARGV != 2;

my $fasta=shift;
my $nlr=shift;



my ($seqID, %seq);
open FA, "$fasta" or die "";
while (<FA>){
	chomp;
	#CNls NLs RNLs TNLs
	if (m/^>(\S+)/){$seqID=$1;}
	else{$seq{$seqID}.=$_;}
}
close FA;

my %nlr;
my @major = ('1', '4', '3', '7');
my @check = ("1437", "143", "147", "137", "437");
open IN, "$nlr" or die "";
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	my @motif=split/,/, $line[-1];

	my $major= 'n';
	my $major_motif="";
        foreach my $motif (@motif){
                foreach my $major (@major){
			$major_motif.="$major"  if ($motif eq $major);
		}
	}
	#print "$major_motif\n";
	next if $major_motif eq "";
	foreach my $check (@check){
		$major = 'y' if $major_motif =~ m/$check/;
	}

	next unless $major eq 'y';


	#my @only = grep{!$motif{$_}} @check;
	#next unless scalar @only == 0;
	$nlr{$line[0]}=1;
	#print "$line[0]\n";
}

close IN;

foreach my $id (keys %seq){
	next if length($seq{$id}) <160;
	my $idd= $id;
	$idd=~ s/\.CNls//g;
        $idd=~ s/\.NLs//g;
        $idd=~ s/\.RNLs//g;
        $idd=~ s/\.TNLs//g;

	next unless exists $nlr{$idd};
	print ">$id\n$seq{$id}\n";
}


