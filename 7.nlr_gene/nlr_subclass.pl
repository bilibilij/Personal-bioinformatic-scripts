#!/usr/bin/envs perl -w

use strict;


my $fasta=shift;
my $tsv=shift;

open FA,"$fasta" or die "";

my ($seqID, %seq);
while (<FA>){
	chomp;
	if (m/^>(\S+)/){$seqID=$1;}
	else{$seq{$seqID}.=$_;}
}
close FA;

my %domain;
open TV, "$tsv" or die "";
while (<TV>){
	chomp;
	my @line=split/\t/, $_;
	next unless (m/RPW8/ || m/PF01582/ || m/Coils/);
	next unless exists $seq{$line[0]};
	next if ($line[7]- $line[6] < 40);
	next if $line[6] >80;
	$domain{$line[0]}{$line[6]}=$line[5];
}
#Coils
close TV;

`rm $fasta.subclass` ;

#open OUT, ">$fasta.subclass" or die "";
foreach my $id (sort keys %seq){
	open OUT, ">>$fasta.subclass" or die "";
	if (not exists $domain{$id} ){print "$id.NLs\n";print OUT  ">$id.NLs\n$seq{$id}\n";}
	else{
		my @pos = (sort {$a<=>$b} keys %{$domain{$id}});
		my $pos=shift @pos;
		if ($domain{$id}{$pos} =~ m/RPW8/){print "$id.RNLs\n";print OUT  ">$id.RNLs\n$seq{$id}\n";}
		elsif ($domain{$id}{$pos} =~ m/TIR/){print "$id.TNLs\n"; print OUT ">$id.TNLs\n$seq{$id}\n";}
		elsif ($domain{$id}{$pos} =~ m/Coil/){print "$id.CNls\n";print OUT ">$id.CNls\n$seq{$id}\n";}
		else{die "there are other domain here$!\n";}
	}
	close OUT;
}
