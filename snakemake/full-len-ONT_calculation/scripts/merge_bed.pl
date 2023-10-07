#!/usr/bin/perl -w

my $dir=shift;
opendir DIR, "$dir" or die "";
my @dir=readdir DIR;

my %bed;
my %check;
foreach my $file (@dir){
	next unless $file =~ m/bed$/;
	open IN, "$dir/$file" or die "";
	while (<IN>){
		chomp;
		my @line =split/\s+/, $_;
		die "err at $dir/$file $!\n" unless @line == 12;
		$bed{$line[0]}{$line[1]}{$line[2]}=$_;
		die "the same ONT reads name between $check{$line[3]} and $file, named  $line[3]\n $! " if exists $check{$line[3]};
		$check{$line[3]}=$file;
	}
	close IN;
}


foreach my $chr (sort {$a cmp $b}  keys %bed){
	foreach my $ss (sort {$a<=>$b} keys %{$bed{$chr}}){
		foreach my $en (sort {$a<=>$b} keys %{$bed{$chr}{$ss}}){
			print "$bed{$chr}{$ss}{$en}\n";
		}
	}
}






