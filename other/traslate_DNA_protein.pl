#!/usr/bin/perl -w

use strict;
use Bio::Seq;

my $usage=<<USAGE;

perl $0 fasta(DNA) > fasta(protein)

USAGE

die "$usage\n" unless @ARGV == 1 ;


#my $fa_in=shift;
open FA, "$ARGV[0]" or die "";

my ($seqID, %seq);
while (<FA>){
	chomp;
	if (m/^>(\S+)/){#$seqID= $1;}
		print "$_\n";
	}
	else{
		my $DNA=$_;
		my $pep=&TranslateDNASeq($DNA);
		#### if you want to translate protein with the "*" replaced by "X", use line of below;
		#$pep =~ s/\*/X/g;
		print "$pep\n";
	}
}

close FA;


#my $pep=&TranslateDNASeq($DNA);

sub TranslateDNASeq{
    use Bio::Seq;
    (my $dna)=@_;
    my $seqobj=Bio::Seq->new(-seq =>$dna, -alphabet =>'dna');
    return $seqobj->translate()->seq();
}






