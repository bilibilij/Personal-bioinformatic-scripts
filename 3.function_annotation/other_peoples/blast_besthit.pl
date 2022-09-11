#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;
usage:
        perl $0 blast_out_fmt6 > blast_out_fmt6.besthit
print best hit by filtering score and identity of blast.tab ;
score > identity : under the setting Evalue

USAGE

if (@ARGV != 1) {die $usage;}

my %stat;
open IN, "$ARGV[0]" or die "can not open file ";
while (<IN>) {
        chomp;
        @_=split /\t/,$_;
                $stat{$_[0]}{$_[11]}{$_[2]}=$_;}
foreach my $pep_id ( keys %stat) {
        my @coverage;
        my @score_blast;
        foreach my $score ( keys %{$stat{$pep_id}} ) {
                push @score_blast, $score;
                @score_blast= sort {$b <=> $a} @score_blast;}
                my $scor = $score_blast[0];
                foreach my $coverage (keys %{$stat{$pep_id}{$scor}} ) {
                                push  @coverage,$coverage;
                                @coverage = sort {$b <=> $a} @coverage;}
                my $cover=$coverage[0];
                print "$stat{$pep_id}{$scor}{$cover}\n";
        }

