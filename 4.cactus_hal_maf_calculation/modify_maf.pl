#!/usr/bin/envs perl -w

use strict;

my $usage =<<USAGE;

perl $0 input.maf > input.maf.modified

USAGE

die "$usage\n" unless @ARGV == 1;




open MF, "$ARGV[0]" or die "";

while (<MF>){
	chomp;
	if (m/^s\s+/){
		my @line=split/\s+/, $_;
		if ($line[1] =~ m/refChr/){
			#print "$_\n";
			next
		}else{
			my @arr=split/\./, $line[1];
			my $last= pop @arr;
			my $pre=join("_", @arr);
			$line[1]="$pre.$last";
			my $print=join("\t", @line);
			print "$print\n";
		}
	}elsif (m/^# hal/){
		print "# hal (c_par:0.365342,(T_hass_prot:6.12269e-18,(A_ara_pep:0.0902432,(Aalp_pep:0.114998,((Mpyg_pep:0.129859,(Lala:0.142895,((C_rub_pep:0.147493,Bstr:0.189519)0.539435:0.0366159,(A_tha_pep:0.113207,Aly_pep:0.0599314)0.794223:0.0110889)0.722006:0.00143956)0.383479:0.00141118)0.148815:1.56741e-16,((T_arv_pep:0.115312,E_sal_pep:0.119472)0.359948:0.0312753,(T_par_pep:0.157187,((S_iri:0.127366,Il_best_pep:0.20511)0.244239:0.00136781,(Ovio_best_pep:0.264257,((Brapa_genome_v3_0_pep:0.556097,BoleraceaHDEM_proteins:0.163903)0.821162:0.0262108,Bnigra_NI100_v2_pep:0.142656)0.851672:0.00130054)0.155144:0.0037011)0.181759:0.0018251)0.289354:0.00126716)0.33236:0.0204311)0.156605:0.0360732)0.709185:0.00126504)0.87926:0.294433)1:0.365342)Anc00;\n";
	}
	
	
	else{
		print "$_\n";
		next;
	}

}


close MF;
