#!/usr/bin/perl  -w

use strict;

my $usage=<<USAGE;

perl $0 meth_bin_intersect bin_num out_file

USAGE

die $usage if @ARGV == 0 ;

#加个第三列， 组的这个参数；


open IN, "$ARGV[0]" or die "";
my $num=$ARGV[1]+1;
my $out=$ARGV[2];



my %bin;
my %number;
while(<IN>){
    chomp;
    my @line=split/\t/,$_;
    my $meth=$line[7];
    my @a=split"_", $line[3];
    if (@a == 4){
	if ($a[1] eq "+"){
		$bin{$a[2]}{$a[3]}+=$meth;
		$number{$a[2]}{$a[3]}++;
	}elsif($a[1] eq "-"){
		my $r=$num-$a[-1];
		$bin{$a[2]}{$r}+=$meth;
		$number{$a[2]}{$r}++;
	}else{die "error \n";}
    }elsif (@a == 3 ){
	if ($a[1] eq "+"){
		$bin{"gene"}{$a[-1]}+=$meth;
		$number{"gene"}{$a[-1]}++;
	}elsif($a[1] eq "-"){
		my $r=$num-$a[-1];
		$bin{"gene"}{$r}+=$meth;
		$number{"gene"}{$r}++;
	}
    }else{die "err\n";}

}

#Chr01   5897    5898    Ovio00001_-_10  Chr01   5897    5898    0       0       3
#Chr01   5979    5980    Ovio00001_-_up_1        Chr01   5979    5980    33.3333333333333        1       3



foreach my $type (keys %bin){
	foreach my $bin (keys %{$bin{$type}}){
		my $c;
		if ($type eq "gene"){
			$c=$ARGV[1]+$bin;
			
		}elsif ($type eq "up"){
			$c=$bin;
		}elsif ($type eq "down"){
			$c=$bin+$ARGV[1]+$ARGV[1];
		}else{die "err";}
		my $d=sprintf "%.2f", ($bin{$type}{$bin}/$number{$type}{$bin});	#$$var=sprintf "%.2f",$var;
	
		print "$c\t$d\t$out\n";
		#print "$c\t$bin{$type}{$bin}\t$out\n";
	}
}




#Chr01   5897    5898    Ovio00001_-_10  Chr01   5897    5898    0       0       3
#Chr01   5979    5980    Ovio00001_-_up_1        Chr01   5979    5980    33.3333333333333        1       3


