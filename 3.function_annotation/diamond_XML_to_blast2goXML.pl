#!/usr/bin/perl 

use strict;

my $usage= <<USAGE;
Usage:
	blast2GOv2.5是老版本用于GO注释的软件，其设计时对应的NR版本也比较低。在2017年之前的nr/nt数据库是支持gi号搜索的，在2017年之后的nr/nt数据库变成不再支持gi号搜索的，变成accession id.
	我们现在做的NR数据库比对，得到的结果都是只有accession id的，因此需要在结果中加入gi号，处理成blast2GO软件能够读取的格式，“prot.accession2taxid.gz”存储着所有accesion id 与gi号的对应关系。
	最好下载最新版本的prot.accession2taxid.gz, 当然还是会有些accessions 没有对应的gi号，不过对于Blast2go软件来说，少一些blast结果影响不大。
	!!!blast2go软件一次性仅支持最多1000条序列进行go注释，需要分开跑!!!
贾畅夫
USAGE


open XML, "$ARGV[0]" or die "can't open diamond XML file";
my %hash=();

while (<XML>){
	chomp;
	if (m/<Hit_id>(\S+)<\/Hit_id>/){
		my $id=$1;
		$hash{$id}=1;
	}
}

close XML;


my $prot=$ARGV[1];
open (IND, "gzip -dc $prot|") or die "can't open prot.accession2taxid.gz"  ;
my %gii=();

while (<IND>){
	chomp;
	my $acc=(split /\t/, $_)[1];
	my $gi=(split /\t/, $_)[3];
	if (exists $hash{$acc}) { 
		$gii{$acc}="gi|$gi|ref|$acc";
	}
}

#gi|118484601|ref|ABK94174.1

close IND;
		


open XML, "$ARGV[0]" or die "can't open diamond XML file";

while (<XML>){
	chomp;
	if (m/BlastOutput_version/){
		print "  <BlastOutput_version>BLASTP 2.2.28</BlastOutput_version>\n";
	}elsif (m/<Hit_id>(\S+)<\/Hit_id>/){
		my $id=$1;
		my $giid=$gii{$id};
		print "\t<Hit_id>$giid</Hit_id>\n";
	}else{
		print "$_\n";
	}
}
