#!/usr/bin/envs perl -w

use strict;

my $usage=<<USAGE;

perl $0 maf

script will create each splited maf by chromosomes of reference at the path of maf input.

USAGE

die "$usage\n" unless @ARGV == 1;


#my @out=split/\./, $ARGV[0];
#pop @out;
#my $out_dir=join("\t", @out);
my $out_dir=$ARGV[0];
$out_dir=~s/\./_/g;
mkdir $out_dir;

open MF, "$ARGV[0]" or die "";
$/="\na\n";
my ($header, %hash);
while (<MF>)
{
        chomp;
        if (m/^#/)
        {
                $header="$_\n";
        }else{
                my @line=split/\s+/, $_;
                $hash{$line[1]}{$line[2]}=$_;
        }
}

close MF;

$/="\n";

foreach my $chr (sort keys %hash){
	#next if ($chr =~ m/ChrM/ || $chr =~m/ChrC/);
	#next unless $chr =~ m/A_tha_pep/;
        my $out=split/\./, $chr;
        open OUT, ">$out_dir/$chr.maf" or die "";
	print OUT "$header\n";
        foreach my $pos (sort {$a<=>$b} keys %{$hash{$chr}}){
                print OUT "a\n$hash{$chr}{$pos}";
        }
        close OUT;
        print "chr has been done! \n ";
}


