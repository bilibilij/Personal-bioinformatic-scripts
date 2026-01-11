#!/usr/bin/perl -w
use strict;

my $usage = <<USAGE;

perl $0 gff(evm or pasa output for sorting) gene_prefix > final.gff

contact jcf with 1020160171\@qq.com

USAGE

die $usage if @ARGV == 0;

open GF, "$ARGV[0]" or die "";

my $h = $ARGV[1]; # 基因前缀
my $m = 0;
my %hash;
my %gene;
my $start;
my $end;
my $gene_mark = 0;
my (@chr, @sca); #store the chr and scaffold sepearately.
while (<GF>) {
    chomp;
    next if m/^$/;
    next if m/^#/;
    my @line = split /\t/, $_;
    pop @line;
    die $! if @line != 8;
    my $c = join("\t", @line);

    if ($line[2] eq "gene") {
        $gene_mark++;
        $start = $line[3];
        $end = $line[4];
        $m = 0;
        $gene{$line[0]}{$start}{$gene_mark} = "$c";
    } elsif ($line[2] eq "mRNA") {
        $m++;
        $hash{$line[0]}{$start}{$gene_mark}{$m} .= "$c\n";
    } else {
        $hash{$line[0]}{$start}{$gene_mark}{$m} .= "$c\n";
    }
}

close GF;

# 自然排序函数
sub natural_sort {
    my @list = @_;
    return sort { 
        my ($anum) = $a =~ /(\d+)/;
        my ($bnum) = $b =~ /(\d+)/;
        $anum <=> $bnum;
    } @list;
}

# 基因 ID 生成规则
#my $gene_num = 100;
my %chromosome_order;
my $chr_index = 1;
my $scaffold_index = 1;

# 先对所有的染色体和 scaffold 进行排序
foreach my $chr (natural_sort(keys %hash)) {
    if ($chr =~ /^scaffold/i) {
        $chromosome_order{$chr} = sprintf "T%03d", $scaffold_index++;
    } else {
        $chromosome_order{$chr} = sprintf "%03d", $chr_index++;
    }
}

###I print chromosome and scaffold repeately and  separately.
#Chr
foreach my $chr (natural_sort(keys %hash)) {
    next if $chr =~ m/scaffold/;
    my $gene_num = 1;
    foreach my $st (sort {$a <=> $b} keys %{$hash{$chr}}) {
        foreach my $ed (sort {$a <=> $b} keys %{$hash{$chr}{$st}}) {
            my $chr_prefix = $chr;
            my $gene_id = sprintf "%s%sG%06d", $h, $chr_prefix, $gene_num;

            # 打印 gene 行
            print "$gene{$chr}{$st}{$ed}\tID=$gene_id;Name=$gene_id;\n";
            $gene_num++;

            # mRNA 和其他特征
            my $mRNA_num = 1;
            foreach my $mr (sort {$a <=> $b} keys %{$hash{$chr}{$st}{$ed}}) {
                my $mRNA = "$gene_id.t$mRNA_num";
                die "err\n" unless my @arr = split /\n/, $hash{$chr}{$st}{$ed}{$mr};
                my %num;

                foreach my $fea (@arr) {
                    my @line = split /\t/, $fea;

                    if ($line[2] eq "mRNA") {
                        print "$fea\tID=$mRNA;Parent=$gene_id;Name=$mRNA;\n";
                    } else {
                        $num{$line[2]}++;
                        my $N = $num{$line[2]};
                        print "$fea\tID=$mRNA.$line[2]$N;Parent=$mRNA;\n";
                    }
                }
                $mRNA_num++;
            }
            print "\n";
        }
    }
}

my $gene_num = 100;
#gene_num = 1;
#scaffold
# 根据排序后的规则生成基因 ID
foreach my $chr (natural_sort(keys %hash)) {
    next unless $chr =~ m/scaffold/;
    foreach my $st (sort {$a <=> $b} keys %{$hash{$chr}}) {
        foreach my $ed (sort {$a <=> $b} keys %{$hash{$chr}{$st}}) {
            my $chr_prefix = $chr;
            my $gene_id = sprintf "%sUn%06d", $h, $gene_num;

            # 打印 gene 行
            print "$gene{$chr}{$st}{$ed}\tID=$gene_id;Name=$gene_id;\n";
            $gene_num++;

            # mRNA 和其他特征
            my $mRNA_num = 1;
            foreach my $mr (sort {$a <=> $b} keys %{$hash{$chr}{$st}{$ed}}) {
                my $mRNA = "$gene_id.t$mRNA_num";
                die "err\n" unless my @arr = split /\n/, $hash{$chr}{$st}{$ed}{$mr};
                my %num;

                foreach my $fea (@arr) {
                    my @line = split /\t/, $fea;

                    if ($line[2] eq "mRNA") {
                        print "$fea\tID=$mRNA;Parent=$gene_id;Name=$mRNA;\n";
                    } else {
                        $num{$line[2]}++;
                        my $N = $num{$line[2]};
                        print "$fea\tID=$mRNA.$line[2]$N;Parent=$mRNA;\n";
                    }
                }
                $mRNA_num++;
            }
            print "\n";
        }
    }
}
