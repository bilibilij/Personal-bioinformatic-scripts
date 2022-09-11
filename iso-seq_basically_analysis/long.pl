#!/usr/bin/perl 
#

open IN, "$ARGV[0]" or die "";
my %seq;
while (<IN>){
	chomp;
	next unless m/^>PB.([^\.]*)\.([^\|]*)\|.*length=(\S+)/;
	my $a=$1;
	my $b=$3;
	my $c=$2;
	my $g="PB.$a";
	my $t="PB.$a.$c";
	$seq{$g}{$t}=$b;
	#print "$a\t$b\n";

}
#>PB.1.1|scaffold0025:225945-231457(+)|transcript/1818 transcript/1818 full_length_coverage=2;length=4064
close IN;

my %best;
foreach my $gene (sort keys %seq){
        my $len=0;
        my $best_trans;
        foreach my $trans ( sort keys %{$seq{$gene}}){
                #for debug, you can wathc this information to check this script is correct or not
                #my $a=length($seq{$gene}{$trans});
                #print ">$gene.$trans\t$a\n";
                #print ">$gene.$trans\n$seq{$gene}{$trans}\n";
                if ($seq{$gene}{$trans} > $len) { $len = $seq{$gene}{$trans}; $best_trans = $trans;}
        }
	#push @best, $best_trans;
	$best{$best_trans}=1;
	#print "$gene\t$best_trans\n";
}

open IN, "$ARGV[1]" or die "";
#my %match;
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	print "$_\n" if exists $best{$line[0]};

}

close IN;






