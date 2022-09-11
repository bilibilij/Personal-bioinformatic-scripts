#!/usr/bin/perl -w 
use strict;

my $cpu=$ARGV[2];
my $mem=$ARGV[3];
my $fasta=`basename $ARGV[0]`;
my @file_name = split /\./, $fasta;
my$out_name=$file_name[0];

`rm -r $out_name.tmp` if (-e "$out_name.tmp");
mkdir "$out_name.tmp" ;

my $mark=1;
my $seq_number=$ARGV[1];
my $num_loop=1;
open OUT, ">$out_name.tmp/$out_name.$mark.fasta" or die "can not open $!";
open IN, "$ARGV[0]" or die "can not open pepFasta";
$/=">"; <IN>; $/="\n";
while (<IN>){
	my $seq_name=">".$_;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$/="\n";
	if ($num_loop % $seq_number ==0 ) {
		$mark++;
		open OUT, ">$out_name.tmp/$out_name.$mark.fasta" or die "can not open $!";}
	print OUT $seq_name.$seq;
	close (OUT) if ( ($num_loop % $seq_number) == ($seq_number-1) );
	$num_loop++;
}


`rm -r $out_name.command.tmp` if (-e "$out_name.command.tmp");
mkdir "$out_name.command.tmp" ;

open COMMAND, ">$out_name.command.tmp/$out_name.command.list";
#blastn -db /export/databases//NT/nt_20201231/nt -query $id -out ${id%.*} -p 1
for (1..$mark) { 
	print COMMAND "/export/software/miniconda3/envs/blast/bin/blastn -db /export/databases//NT/nt_202006/nt  -query ../$out_name.tmp/$out_name.$_.fasta  -out ../$out_name.tmp/$out_name.$_.ntBlast.tab -num_threads  $cpu -outfmt 5  -evalue 1e-5 -max_target_seqs 20 \n";
}
chdir "$out_name.command.tmp";
my $job_mark=1;
open JOB, "$out_name.command.list";
while (<JOB>){
	chomp;
	open SUBMIT, ">R$out_name.$job_mark.sh";
	print SUBMIT "$_\n";
	system `qsub -pe smp $cpu -l vf=$mem R$out_name.$job_mark.sh`;
	$job_mark++;
	close SUBMIT;
}

