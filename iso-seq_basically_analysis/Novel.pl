#!/usr/bin/perl -w

use strict;

open IN, "$ARGV[0]" or die "";
my %hash1;
my %trans1;
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	die $!  unless m/PB\.([^.]*)\.(\S+)\t/;
	my $trans = shift @line;
	unshift @line, "PB.$1";
	my $p=join("\t", @line);
	$hash1{$line[0]}=$p;
	$trans1{$line[0]}=$trans;
	#print "$p\n";
}

close IN;


open IN, "$ARGV[1]" or die "";
my %hash2;
my %trans2;
while (<IN>){
        chomp;
        my @line=split/\t/, $_;
        die $! unless m/PB\.([^.]*)\.(\S+)\t/;
        my $trans = shift @line;
        unshift @line, "PB.$1";
        my $p=join("\t", @line);
        $hash2{$line[0]}=$p;
	$trans2{$line[0]}=$trans;

}

close IN;


open IN, "$ARGV[2]" or die "";
my %seq;
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	die $!  unless $line[3] =~ m/PB\.([^.]*)\.(\S+)/;
	my $h1=$1;
        die $!  unless $line[7] =~ m/PB\.([^.]*)\.(\S+)/;
	my $h2=$1;
	$seq{$h1}{$h2} = $line[-1];
}
close IN;


my (%h1, %h2);
foreach my $gene (sort keys %seq){
        my $len=0;
        my $best_trans;
        foreach my $trans ( sort keys %{$seq{$gene}}){
                if ($seq{$gene}{$trans} > $len) { $len = $seq{$gene}{$trans}; $best_trans = $trans;}
        }
	$h1{$gene}=$best_trans;
	$h2{$best_trans}=$gene;
}




my ($seqID, %seq1, %seq2);
open IN, "$ARGV[3]" or die "";
while (<IN>){
	chomp;
	if (m/^>([^|]*)\|/){$seqID=$1;}
	else{$seq1{$seqID}.=$_;}

}
close IN;



open IN, "$ARGV[4]" or die "";
while (<IN>){
        chomp;
        if (m/^>([^|]*)\|/){$seqID=$1;}
        else{$seq2{$seqID}.=$_;}

}
close IN;





my $Novel = 1;
#%trans1, %trans2;
open OUT, ">Novel.fasta" or die "";
foreach my $hg1 (sort keys %hash1){
	#my $value = $hash1{$hg1};
	 if (exists $h1{$hg1} ){
		if (exists  $hash2{$h1{$hg1}}){
			my $id= sprintf "%#07s", $Novel;
        	        $id = "Novelgene$id";
	                $Novel++;
			print "$hash1{$hg1}\tshared\th1\t$id\n";
		}
	 }else{
	         my $id= sprintf "%#07s", $Novel;
                 $id = "Novelgene$id";
                 $Novel++;
		 print "$hash1{$hg1}\tunique\th1\t$id\n";

		 my $ID = $trans1{$hg1};
		 print OUT ">$id|$ID|L596\n$seq1{$ID}\n";
	}
}


foreach my $hg2 (sort keys %hash2){
        #my $value = $hash1{$hg1};
         if (exists $h2{$hg2} ){
                if (exists  $hash1{$h2{$hg2}}){
                        my $id= sprintf "%#07s", $Novel;
                        $id = "Novelgene$id";
                        $Novel++;

                        print "$hash2{$hg2}\tshared\th2\t$id\n";
		}

         }else{
                 my $id= sprintf "%#07s", $Novel;
                 $id = "Novelgene$id";
                 $Novel++;
                 print "$hash2{$hg2}\tunique\th2\t$id\n";

                 my $ID = $trans2{$hg2};
                 print OUT ">$id|$ID|H596\n$seq2{$ID}\n";


        }
}


close OUT;


