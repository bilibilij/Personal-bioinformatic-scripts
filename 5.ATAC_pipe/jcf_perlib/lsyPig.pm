#!/usr/bin/env perl -w

use strict;
use Exporter;
package lsyPig;

BEGIN{
        print "=============Hello everybody, I'm LSY, a good \"?\" from SiChuan Univiersity, my research
                \"Demographic History and Natural Selection Shape Patterns of Deleterious Mutation Load and Barriers to Introgression across Populus Genome\"
                has been published on Molecular Biology and Evolution(Feburary 03,2022),and I like eating shit, welcome to join me!~ ===================\n\n\n";

}END{
        print "\n\n\n=============If you are interseted in population TE, SV, introgression, GWAS, NGS population analysis,sampling collection, drinking coffee,
       	please contact me, LSY, YYDS, 13161998855\@qq.com or my pig WJL 1053364212\@qq.com===================\n";

}


sub read_fa_from_dir{
	my $seqID=shift;
	my $seq=shift;
	my $files=shift;
	my $dir=shift;
	print "====== read fasta beginning ========\n";
	foreach my $file (@$files){
	        next unless $file =~ m/fa$/;
	        open FA, "$dir/$file" or die "";
	        while(<FA>){
	                chomp;
	                if (m/^>(\S+)/){$seqID=$1;}
	                else{$$seq{$seqID}.=$_;}
	        }
	        close FA;
	}
	print "====== read fasta finished ========\n";
}


sub read_gene_from_dir{
        my $seq=shift;
        my $files=shift;
        my $dir=shift;
	my $seqID;
        print "====== read fasta beginning ========\n";
        foreach my $file (@$files){
                next unless $file =~ m/fa$/;
                open FA, "$dir/$file" or die "";
                while(<FA>){
                        chomp;
                        if (m/^>(\S+)/){$seqID=$1; my @id=split/\./,$seqID; pop @id; $seqID = join(".",@id);  }
                        else{$$seq{$seqID}.=$_;}
                }
                close FA;
        }
        print "====== read fasta finished ========\n";
}


sub check_dependency {
        my $software_info =shift;
        my $match = shift;
        if ($software_info =~ m/$match/) {
            print STDERR "$match:\tOK\n";
        }
        else {
            die "$match:\tFailed\n\n";
        }
}


sub command {
	my $command=shift;
        print STDERR  (localtime) . ": CMD: $command\n";
        system("$command") == 0 or die  "failed to execute: $command\n";
}

sub get_fa_len {

	my $fa=shift;
	my $match=shift;
	my ($seqID, %seq);

	open IN, "$fa" or die "";
	while (<IN>){
	        if (m/^>(\S+)/){$seqID = $1;}
		else{$seq{$seqID}.=$_;}
	}
	close IN;

	my $len;
	foreach my $id (keys %seq){
	        next unless $id =~m/^$match/;
		#print OUT ">$id\n$seq{$id}\n";
	        $len += length($seq{$id}) ;
	}
	
	return($len);

}



#for directly use function name of this pm;
our @ISA = qw(Exporter);
our @EXPORT = qw/check_dependency read_fa_from_dir read_gene_from_dir command get_fa_len/;



1
