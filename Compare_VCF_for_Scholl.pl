#! perl -w

## This script needs Fisher test module

use lib '/Library/Perl/5.8.8/modules';
use lib '/panfs/home/mc854/programs/Text-NSP-1.23';
use Text::NSP::Measures::2D::Fisher2::twotailed;

my $flag_dbSNP = "/home/mc854/anno_library_hg19/snp135Flagged_hg19.txt";
my $in_house_dir = "/home/mc854/gene_burden/in_house";

my ($tgn_denom, $nhlbi_ex_denom, $yale_ex_denom) = (1094*2, 6500*2, 2000*2);	# for calculating_freq subroutine
my $freq_denominator = $tgn_denom + $nhlbi_ex_denom + $yale_ex_denom;

if (@ARGV == 8) {
	my ($include_db_tg, $qs_cutoff, $freq, $nih, $swed, $ct_vcf, $ca_vcf, $outfile) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4], $ARGV[5], $ARGV[6], $ARGV[7]);
	die "## Need number for QS cutoff!\n" if ($qs_cutoff =~ /[A-Z]/);
	die "## Need number for frequency!\n" if ($freq =~ /[A-Z]/);
	die "## Wrong argument for dbSNP / 1000 genomes!\n" if ($include_db_tg !~ /[yn]/i);
	die "## Wrong argument for NHLBI!\n" if ($nih !~ /[e0-9]/i);
	die "## Wrong argument for in_house controls!\n" if ($swed !~ /[ie]/i);
	$ct_vcf = "${in_house_dir}/in_house_controls.hg19.vcf_anno.vcf" if ($ct_vcf =~ /^s$/i);
	$nih = 0 if ($nih =~ /e/i);

	my ($snv_qscutoff, $indel_qscutoff) = ($qs_cutoff, 10*$qs_cutoff);
	if ($include_db_tg =~ /y/i) {
		print "## Will include dbSNP or 1000 genomes variants, just filter by frequency.\n";
	} else {
		print "## Will NOT include dbSNP or 1000 genomes variants (traditional approach).\n";
	}
	print "## SNV QS cutoff = $snv_qscutoff, Indel QS cutoff = $indel_qscutoff.\n";

	my ($ct_no, $ca_no) = (0, 0);

	my (@table, @ct_variants, @ca_variants);
	my $header = "\t\t\t\t\t\tTest 1. No of missense and damaging alleles\t\t\t\t\tTest 2. No of conserved missense and damaging alleles\t\t\t\t\tTest 3. No of damaging alleles\t\t\t\t\tTest 4. No of HOM missense and damaging alleles\t\t\t\t\tTest 5. No of HOM conserved missense and damaging alleles\t\t\t\t\tTest 6. No of HOM damaging alleles
Gene\tFull name\tOMIM\tSegdup\tNo loci in controls\tNo loci in cases\tA\tB\tC\tD\tp-value\tA\tB\tC\tD\tp-value\tA\tB\tC\tD\tp-value\tA\tB\tC\tD\tp-value\tA\tB\tC\tD\tp-value\tA\tB\tC\tD\tp-value\n";

	push @table, "## Control vcf: $ct_vcf\n## Case vcf: $ca_vcf\n";
	if ($freq == 0) {
		print "## Filter for novel-novel..\n";
		push @table, "## Filter for novel-novel..\n";
	} else {
		print "## Filter by frequency ${freq}% across 1000 genomes, NHLBI exomes and Yale DB..\n";
		push @table, "## Filter by frequency ${freq}% across 1000 genomes, NHLBI exomes and Yale DB..\n";
	}

	my %swed;
	if ($swed =~ /i/i) {
		print "## Filter against in_house high BP group..\n";
		push @table, "## Filter against in_house high BP group..\n";
		open SW, "${in_house_dir}/in_house_cases.hg19.vcf_anno.vcf" or die "## Cannot find in_house high BP cohort to be used as controls!\n";
		while (<SW>) {
			next if ($_ =~ /^#/);
			next if ($_ =~ /^chr(23|25|M|Y)\t/);
			next if ($_ =~ /\t(Intergenic|Intron|[35]UTR|Coding\-silent|Indel)\t/);

			if ($_ =~ /^([^\t]+\t\d+)\t/) {	$swed{$1} = 1;	}
		}
		close SW;
	}

	my ($ct_header, $ca_header);
	open CT, $ct_vcf or die "## Cannot find $ct_vcf!\n";
	while (<CT>) {
		last if ($ct_no > 0);
		if ($_ =~ /^\#CHROM\tPOS\t.+\tFORMAT\t(.+)\tGene\tFull name/) {
			my @c = split /\t/, $1;
			$ct_no = @c;
			$ct_header = "Control\t$_";
		}
	}
	close CT;
	open CA, $ca_vcf or die "## Cannot find $ca_vcf!\n";
	while (<CA>) {
		last if ($ca_no > 0);
		if ($_ =~ /^\#CHROM\tPOS\t.+\tFORMAT\t(.+)\tGene\tFull name/) {
			my @c = split /\t/, $1;
			$ca_no = @c;
			$ca_header = "Case\t$_";
		}
	}
	close CA;
	print "## Detected $ct_no controls and $ca_no cases..\n\n";
	push @table, "## Detected $ct_no controls and $ca_no cases..\n\n";
	push @ct_variants, "## Detected $ct_no controls and $ca_no cases..\n\n";
	push @ca_variants, "## Detected $ct_no controls and $ca_no cases..\n\n";
	push @ct_variants, $ct_header;
	push @ca_variants, $ca_header;

	my %flg_dbsnp;
	my $flg_dbsnp = 0;
	open DB, $flag_dbSNP or die "## Cannot find $flag_dbSNP!\n";
	while (<DB>) {
		$flg_dbsnp{$1} = 1 if ($_ =~ /\t(rs\d+)\t/);
		$flg_dbsnp++;
	}
	close DB;
	print "## Will include $flg_dbsnp clinically associated dbSNPs..\n\n";
	push @table, "## Will include $flg_dbsnp clinically associated dbSNPs..\n";

	push @table, $header;

	my (%ct_protalt, %ct_protaltCONS, %ct_dam, %ct_protalt_HOM, %ct_protaltCONS_HOM, %ct_dam_HOM, %ct_comphet, %ct_comphetCONS, %ct_comphet_dam, %ct_loci, %ct_segdup);
	my (%ca_protalt, %ca_protaltCONS, %ca_dam, %ca_protalt_HOM, %ca_protaltCONS_HOM, %ca_dam_HOM, %ca_comphet, %ca_comphetCONS, %ca_comphet_dam, %ca_loci, %ca_segdup);
	my %all_protalt;
	my (%gene_name, %omim);

	my ($ct_snv, $ct_indel) = (0, 0);
	my ($ca_snv, $ca_indel) = (0, 0);

	print "\tReading control variants..\n";
	open CT, $ct_vcf or die "## Cannot find $ct_vcf!\n";
	while (<CT>) {
		next if ($_ =~ /^#/);
		next if ($_ =~ /^chr(23|25|M|Y)\t/);
		next if ($_ =~ /\t(Intergenic|Intron|[35]UTR|Coding\-silent|Indel)\t/);

		my @c = split /\t/, $_;
		my ($position, $ref, $alt, $qs, $gene, $gene_name, $dbsnp, $nhlbi, $tgenome, $yale_freq, $type, $AAcons, $omim, $indel_eval) = ("$c[0]\t$c[1]", $c[3], $c[4], $c[5], $c[8 + $ct_no + 1], $c[8 + $ct_no + 2], $c[8 + $ct_no + 7], $c[8 + $ct_no + 8], $c[8 + $ct_no + 9], $c[8 + $ct_no + 10], $c[8 + $ct_no + 12], $c[8 + $ct_no + 17], $c[8 + $ct_no + 23], indel_eval($c[3], $c[4]));
		$AAcons = 0 if ($AAcons eq "na");
		next if ($qs eq ".");
		if ($include_db_tg =~ /n/i) {
			next if ($dbsnp ne "Novel" && !$flg_dbsnp{$dbsnp});
			next if ($tgenome ne "Novel");
		}
		if ($nhlbi ne "Novel") {
			next if ($nih < $nhlbi);
		}	# instead of ## next if ($nih =~ /i/i && $nhlbi ne "Novel");
		next if ($swed{$position});

		if (length($ref) * length($alt) == 1) {			#snv
			next if ($qs < $snv_qscutoff);
			next if (calculating_freq($tgenome, $nhlbi, $yale_freq) > $freq/100);
			$ct_snv++;
		} else {						#indel
			next if ($qs < $indel_qscutoff);
			next if ($yale_freq > $freq/100);
			$ct_indel++;
		}

		$gene_name{$gene} = $gene_name;
		$omim{$gene} = $omim;
		my ($allele_no, $het_ind, $hom_ind) = (0, 0, 0);
		for my $i (8+1..8+$ct_no) {
			my $var = 0;
			if ($c[$i] =~ /^0[\|\/]1/ || $c[$i] =~ /^1[\|\/]0/) {
				$allele_no++;
				$het_ind++;
				$var = 1;
			} elsif ($c[$i] =~ /^1[\|\/]1/) {
				$allele_no += 2;
				$hom_ind++;
				$var = 1;
			}
			if ($var == 1) {
				$ct_comphet{$gene}{$i}++;
				$ct_comphetCONS{$gene}{$i}++ if ($AAcons <= 1 || $type ne "Coding-missense");
				$ct_comphet_dam{$gene}{$i}++ if ($type eq "Coding-nonsense" || $type eq "Ex-In boundary" || $indel_eval eq "D");
			}
		}

		$ct_loci{$gene}++;
		if ($_ !~ /\tNone\tNone\t/) {	$ct_segdup{$gene} = "Segdup";	}
		$ct_protalt{$gene} += $allele_no;
		$all_protalt{$gene} += $allele_no;
		$ct_protalt_HOM{$gene} += $hom_ind * 2;
		if ($AAcons <= 1 || $type ne "Coding-missense") {
			$ct_protaltCONS{$gene} += $allele_no;
			$ct_protaltCONS_HOM{$gene} += $hom_ind * 2;
		}
		if ($type eq "Coding-nonsense" || $type eq "Ex-In boundary" || $indel_eval eq "D") {
			$ct_dam{$gene} += $allele_no;
			$ct_dam_HOM{$gene} += $hom_ind * 2;
		}

		push @ct_variants, "Control\t$_";
	}
	close CT;

	print "\tReading case variants..\n\n";
	open CA, $ca_vcf or die "## Cannot find $ca_vcf!\n";
	while (<CA>) {
		next if ($_ =~ /^#/);
		next if ($_ =~ /^chr(23|25|M|Y)\t/);
		next if ($_ =~ /\t(Intergenic|Intron|[35]UTR|Coding\-silent|Indel)\t/);

		my @c = split /\t/, $_;
		my ($position, $ref, $alt, $qs, $gene, $gene_name, $dbsnp, $nhlbi, $tgenome, $yale_freq, $type, $AAcons, $indel_eval) = ("$c[0]\t$c[1]", $c[3], $c[4], $c[5], $c[8 + $ca_no + 1], $c[8 + $ca_no + 2], $c[8 + $ca_no + 7], $c[8 + $ca_no + 8], $c[8 + $ca_no + 9], $c[8 + $ca_no + 10], $c[8 + $ca_no + 12], $c[8 + $ca_no + 17], indel_eval($c[3], $c[4]));
		$AAcons = 0 if ($AAcons eq "na");
		next if ($qs eq ".");
		if ($include_db_tg =~ /n/i) {
			next if ($dbsnp ne "Novel" && !$flg_dbsnp{$dbsnp});
			next if ($tgenome ne "Novel");
		}
#		next if ($yale_freq > $freq/100);
		if ($nhlbi ne "Novel") {
			next if ($nih < $nhlbi);
		}	# instead of ## next if ($nih =~ /i/i && $nhlbi ne "Novel");
		next if ($swed{$position});

		if (length($ref) * length($alt) == 1) {			#snv
			next if ($qs < $snv_qscutoff);
			next if (calculating_freq($tgenome, $nhlbi, $yale_freq) > $freq/100);
			$ca_snv++;
		} else {						#indel
			next if ($qs < $indel_qscutoff);
			next if ($yale_freq > $freq/100);
			$ca_indel++;
		}

		$gene_name{$gene} = $gene_name;
		my ($allele_no, $het_ind, $hom_ind) = (0, 0, 0);
		for my $i (8+1..8+$ca_no) {
			my $var = 0;
			if ($c[$i] =~ /^0[\|\/]1/ || $c[$i] =~ /^1[\|\/]0/) {
				$allele_no++;
				$het_ind++;
				$var = 1;
			} elsif ($c[$i] =~ /^1[\|\/]1/) {
				$allele_no += 2;
				$hom_ind++;
				$var = 1;
			}
			if ($var == 1) {
				$ca_comphet{$gene}{$i}++;
				$ca_comphetCONS{$gene}{$i}++ if ($AAcons <= 1 || $type ne "Coding-missense");
				$ca_comphet_dam{$gene}{$i}++ if ($type eq "Coding-nonsense" || $type eq "Ex-In boundary" || $indel_eval eq "D");
			}
		}

		$ca_loci{$gene}++;
		if ($_ !~ /\tNone\tNone\t/) {	$ca_segdup{$gene} = "Segdup";	}
		$ca_protalt{$gene} += $allele_no;
		$all_protalt{$gene} += $allele_no;
		$ca_protalt_HOM{$gene} += $hom_ind * 2;
		if ($AAcons <= 1 || $type ne "Coding-missense") {
			$ca_protaltCONS{$gene} += $allele_no;
			$ca_protaltCONS_HOM{$gene} += $hom_ind * 2;
		}
		if ($type eq "Coding-nonsense" || $type eq "Ex-In boundary" || $indel_eval eq "D") {
			$ca_dam{$gene} += $allele_no;
			$ca_dam_HOM{$gene} += $hom_ind * 2;
		}

		push @ca_variants, "Case\t$_";
	}
	close CT;

	foreach my $gene (keys %ct_comphet) {
		foreach my $sample (keys %{$ct_comphet{$gene}}) {
			$ct_comphetCONS{$gene}{$sample} = 0 if (!$ct_comphetCONS{$gene}{$sample});
			$ct_comphet_dam{$gene}{$sample} = 0 if (!$ct_comphet_dam{$gene}{$sample});
			if ($ct_comphet{$gene}{$sample} >= 2) {	$ct_protalt_HOM{$gene} += 2;	}
			if ($ct_comphetCONS{$gene}{$sample} >= 2) {	$ct_protaltCONS_HOM{$gene} += 2;	}
			if ($ct_comphet_dam{$gene}{$sample} >= 2) {	$ct_dam_HOM{$gene} += 2;	}
		}
	}
	foreach my $gene (keys %ca_comphet) {
		foreach my $sample (keys %{$ca_comphet{$gene}}) {
			$ca_comphetCONS{$gene}{$sample} = 0 if (!$ca_comphetCONS{$gene}{$sample});
			$ca_comphet_dam{$gene}{$sample} = 0 if (!$ca_comphet_dam{$gene}{$sample});
			if ($ca_comphet{$gene}{$sample} >= 2) {	$ca_protalt_HOM{$gene} += 2;	}
			if ($ca_comphetCONS{$gene}{$sample} >= 2) {	$ca_protaltCONS_HOM{$gene} += 2;	}
			if ($ca_comphet_dam{$gene}{$sample} >= 2) {	$ca_dam_HOM{$gene} += 2;	}
		}
	}

	foreach my $k (keys %all_protalt) {
		$ct_segdup{$k} = "None" if (!$ct_segdup{$k});
		$ct_loci{$k} = "0" if (!$ct_loci{$k});
		$ct_protalt{$k} = 0 if (!$ct_protalt{$k});
		$ct_protaltCONS{$k} = 0 if (!$ct_protaltCONS{$k});
		$ct_dam{$k} = 0 if (!$ct_dam{$k});
		$ct_protalt_HOM{$k} = 0 if (!$ct_protalt_HOM{$k});
		$ct_protaltCONS_HOM{$k} = 0 if (!$ct_protaltCONS_HOM{$k});
		$ct_dam_HOM{$k} = 0 if (!$ct_dam_HOM{$k});

		$ca_segdup{$k} = "None" if (!$ca_segdup{$k});
		$ca_loci{$k} = "0" if (!$ca_loci{$k});
		$ca_protalt{$k} = 0 if (!$ca_protalt{$k});
		$ca_protaltCONS{$k} = 0 if (!$ca_protaltCONS{$k});
		$ca_dam{$k} = 0 if (!$ca_dam{$k});
		$ca_protalt_HOM{$k} = 0 if (!$ca_protalt_HOM{$k});
		$ca_protaltCONS_HOM{$k} = 0 if (!$ca_protaltCONS_HOM{$k});
		$ca_dam_HOM{$k} = 0 if (!$ca_dam_HOM{$k});

		my ($fish1, $fish2, $fish3, $fish4, $fish5, $fish6) = (1, 1, 1, 1, 1, 1);
		$fish1 = fish($ct_protalt{$k}, 2*$ct_no, $ct_protalt{$k} + $ca_protalt{$k}, 2*$ct_no + 2*$ca_no);
		$fish2 = fish($ct_protaltCONS{$k}, 2*$ct_no, $ct_protaltCONS{$k} + $ca_protaltCONS{$k}, 2*$ct_no + 2*$ca_no);
		$fish3 = fish($ct_dam{$k}, 2*$ct_no, $ct_dam{$k} + $ca_dam{$k}, 2*$ct_no + 2*$ca_no);
		$fish4 = fish($ct_protalt_HOM{$k}, 2*$ct_no, $ct_protalt_HOM{$k} + $ca_protalt_HOM{$k}, 2*$ct_no + 2*$ca_no);
		$fish5 = fish($ct_protaltCONS_HOM{$k}, 2*$ct_no, $ct_protaltCONS_HOM{$k} + $ca_protaltCONS_HOM{$k}, 2*$ct_no + 2*$ca_no);
		$fish6 = fish($ct_dam_HOM{$k}, 2*$ct_no, $ct_dam_HOM{$k} + $ca_dam_HOM{$k}, 2*$ct_no + 2*$ca_no);

		$omim{$k} = "None" if (!$omim{$k});
		my $out = sprintf "%s\t%s\t%s\t%s\_%s\t%d\t%d\t%d\t%d\t%d\t%d\t%E\t%d\t%d\t%d\t%d\t%E\t%d\t%d\t%d\t%d\t%E\t%d\t%d\t%d\t%d\t%E\t%d\t%d\t%d\t%d\t%E\t%d\t%d\t%d\t%d\t%E\n",
$k, $gene_name{$k}, $omim{$k}, $ct_segdup{$k}, $ca_segdup{$k},
$ct_loci{$k}, $ca_loci{$k},
$ct_protalt{$k}, 2*$ct_no - $ct_protalt{$k}, $ca_protalt{$k}, 2*$ca_no - $ca_protalt{$k}, $fish1,
$ct_protaltCONS{$k}, 2*$ct_no - $ct_protaltCONS{$k}, $ca_protaltCONS{$k}, 2*$ca_no - $ca_protaltCONS{$k}, $fish2,
$ct_dam{$k}, 2*$ct_no - $ct_dam{$k}, $ca_dam{$k}, 2*$ca_no - $ca_dam{$k}, $fish3,
$ct_protalt_HOM{$k}, 2*$ct_no - $ct_protalt_HOM{$k}, $ca_protalt_HOM{$k}, 2*$ca_no - $ca_protalt_HOM{$k}, $fish4,
$ct_protaltCONS_HOM{$k}, 2*$ct_no - $ct_protaltCONS_HOM{$k}, $ca_protaltCONS_HOM{$k}, 2*$ca_no - $ca_protaltCONS_HOM{$k}, $fish5,
$ct_dam_HOM{$k}, 2*$ct_no - $ct_dam_HOM{$k}, $ca_dam_HOM{$k}, 2*$ca_no - $ca_dam_HOM{$k}, $fish6;
		push @table, $out;
	}

	print "\tFrom controls, $ct_snv SNVs and $ct_indel indels..\n";
	print "\tFrom cases, $ca_snv SNVs and $ca_indel indels..\n";
	push @ct_variants, "\tFrom controls, $ct_snv SNVs and $ct_indel indels..\n";
	push @ca_variants, "\tFrom cases, $ca_snv SNVs and $ca_indel indels..\n";

	print "## Writing ${outfile}.table.txt..\n";
	open OUT, ">${outfile}.table.txt";
	print OUT @table;
	close OUT;

	print "## Writing ${outfile}.ct_variants.txt..\n";
	open OUT, ">${outfile}.ct_variants.txt";
	print OUT @ct_variants;
	close OUT;
	print "## Writing ${outfile}.ca_variants.txt..\n";
	open OUT, ">${outfile}.ca_variants.txt";
	print OUT @ca_variants;
	close OUT;
} else {
	die "## perl p02Compare_VCF <Include_dbSNP_tgVar_[Y]es_or_[N]o> <qs_cutoff> <freq:[0]_or_[1]%_or_any%> <NHLBI:[E]xclude_or_[Freq]&> <in_house_HighBP:[I]nclude_or_[E]xclude> <control_anno.VCF_or_use_[S]wedish_LowBP_control> <case_anno.VCF> <output_without_extension>\n";
}

sub indel_eval {
	my ($ref, $alt) = ($_[0], $_[1]);
	if ($alt =~ /[ID](\d+)$/) {
		if ($1 % 3 == 0) {
			return "C";
		} else {
			return "D";
		}
	} elsif ($alt =~ /[ID](\w+)$/) {
		if (length($1) % 3 == 0) {
			return "C";
		} else {
			return "D";
		}
	} elsif (length($ref) * length($alt) != 1) {			#indel
		if (abs(length($ref) - length($alt)) % 3 == 0) {			#inframe
			return "C";
		} else {									#outofframe
			return "D";
		}
	} else {
		return "SNV";
	}
}
sub fish {	# ($cell11, $cell1t, $cellt1, $celltt)
	return calculateStatistic(n11=>$_[0], n1p=>$_[1], np1=>$_[2], npp=>$_[3]);
}
sub calculating_freq {
	my $numerator = 0;
	my ($tgn, $nhlbi_ex, $yale_ex) = ($_[0], $_[1], $_[2]);
	if ($tgn =~ /1\:(\d+)\,2\:(\d+)/) {	$numerator += $1 + 2*$2;	}
	if ($nhlbi_ex =~ /\.\d+/) {	$numerator += int($nhlbi_ex * $nhlbi_ex_denom/100);	}
	if ($yale_ex =~ /\.\d+/) {	$numerator += int($yale_ex * $yale_ex_denom);	}
	return $numerator/$freq_denominator;
}
