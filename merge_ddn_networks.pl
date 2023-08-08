#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

sub read_names($$);
sub read_sif_files($$);
sub read_annotations($$);
sub merge_annotations($$);

my $names_array;
my $sifs_array;
my $annotations_array;
my $enrichments_array;
my $thres;
my $out_file;

GetOptions (
		#"help"				=> \$help,			#TO DO	
		"name=s@"			=> \$names_array,
		"sif=s@"			=> \$sifs_array,
		"enrich=s@" 		=> \$enrichments_array,
		"anno=s@" 			=> \$annotations_array,
		"thres=s"			=> \$thres,
		"out=s"				=> \$out_file,
		);

my $data;
$data = read_names($data, $names_array);

$data = read_sif_files($data, $sifs_array);

$data = read_annotations($data, $annotations_array);

my $annotation_output = merge_annotations($data, $thres);
open(OUT , ">", $out_file . '_anno.tsv') or die;
print OUT $annotation_output;
close(OUT);

open(OUT , ">", $out_file . '.sif') or die;
#print OUT "PACA\tINTERACTION\tPACB\n";
foreach my $interaction (sort keys %{$data->{'sif'}{'all'}}) {
	print OUT $interaction . "\n";
}
close(OUT);

########################################################
sub read_names($$) {
	my ($data, $names_array) = @_;
	foreach my $name (@$names_array) {
		push(@{$data->{'name'}}, $name);
	}
	return $data;
}


sub read_sif_files($$) {
	my ($data, $sifs_array) = @_;
	for(my $i = 0; $i < scalar(@$sifs_array); $i++) {
		my $counter = 0;
		open(IN, $sifs_array->[$i]) or die;
		my $name = $data->{'name'}[$i];
		my $header;
		while(defined(my $line = <IN>)) {
			# unless($header) {
			# 	$header++;
			# 	next;
			# }
			chomp($line);
			$counter++;
			$data->{'sif'}{'all'}{$line}++;
			$data->{'sif'}{$name}{$line}++;
			my ($a, $inter, $b) = split(/\t/, $line);
			$data->{'pac'}{$name}{$a}++;
			$data->{'pac'}{$name}{$b}++;
			$data->{'pac'}{'all'}{$a}++;
			$data->{'pac'}{'all'}{$b}++;
		}
		print STDERR $counter . " entries in " . $sifs_array->[$i] . "\n";
		close(IN);
	}
	print STDERR scalar(keys %{$data->{'sif'}{'all'}}) . " entries combined\n";
	return $data;
}

sub read_annotations($$) {
	my ($data, $annotations_array) = @_;
	for(my $i = 0; $i < scalar (@$annotations_array); $i++) {
		my $name = $data->{'name'}[$i];
		open(IN, $annotations_array->[$i]) or die;
		my $header;
		#	0	PAC	
		#	1	gene	
		#	2	bozer_score	
		#	3	mutsigcv_qvalue	
		#	4	hotspots_min_qvalue	
		#	5	hotspots_num_hotspots	
		#	6	study_size	
		#	7	bozer_significant	
		#	8	mutsigcv_significant	
		#	9	hotspots_significant	
		#	10	3of3_significant	
		#	11	2of3_significant	
		#	12	cohort	
		#	13	num_patients_with_nonsyn	
		#	14	num_patients_with_null	
		#	15	num_patients_with_gain	
		#	16	num_patients_with_loss	
		#	17	num_patients_with_snv	
		#	18	ratio_patients_with_snv	
		#	19	num_patients_with_cnv	
		#	20	ratio_patients_with_cnv	
		#	21	num_patients_with_either	
		#	22	ratio_patients_with_either	
		#	23	2of3 and 1pc	
		#	24	2of3 and 5pc	
		#	25	2of3 and 10pc	
		#	26	2of3 and 15pc	
		#	27	2patients	
		#	28	3patients
		#   29  hitlist
		
		while(defined(my $line = <IN>)) {
			unless($header) {
				$header++;
				next;
			}
			chomp($line);
			my @cols = split(/\t/, $line);
			$data->{'anno'}{$cols[0]}{'gene'} = $cols[1];
			$data->{'anno'}{$cols[0]}{'ratio_patients_with_snv'}{$name} = $cols[18];
			$data->{'anno'}{$cols[0]}{'ratio_patients_with_cnv'}{$name} = $cols[20];
			$data->{'anno'}{$cols[0]}{'ratio_patients_with_either'}{$name} = $cols[22];
			$data->{'anno'}{$cols[0]}{'cohort'}{$name} = $cols[12];
			$data->{'anno'}{$cols[0]}{'2of3'}{1}{$name} = $cols[23];
			$data->{'anno'}{$cols[0]}{'2of3'}{5}{$name} = $cols[24];
			$data->{'anno'}{$cols[0]}{'2of3'}{10}{$name} = $cols[25];
			$data->{'anno'}{$cols[0]}{'2of3'}{15}{$name} = $cols[26];
			$data->{'anno'}{$cols[0]}{'hitlist'}{$name} = $cols[29];
			$data->{'anno'}{$cols[0]}{'typeofhit'}{$name} = $cols[33];
		}
		close(IN);
	}
	return $data;
}

sub merge_annotations($$) {
	my ($data, $thres) = @_;
	my @pacs = sort keys %{$data->{'pac'}{'all'}};
	my $output = "PAC\tgene";
	foreach my $name (@{$data->{'name'}}) {
		$output .= "\tcohort " . $name;
		$output .= "\t2of3hit " . $name;
		$output .= "\thitlist " . $name;
		$output .= "\ttypeofhit " . $name;
		$output .= "\tratio_patients_with_either " . $name;
		$output .= "\tvisual " . $name;
		$output .= "\tratio_patients_with_snv " . $name;
		$output .= "\tvisual_snv " . $name;
		$output .= "\tratio_patients_with_cnv " . $name;
		$output .= "\tvisual_cnv " . $name;
	}
	$output .= "\n";
	foreach my $pac (@pacs) {
		my $gene = $data->{'anno'}{$pac}{'gene'};
		unless($gene) {
			$gene = 'FIXME';
		}
		$output .= $pac . "\t" . $gene;

		foreach my $name (@{$data->{'name'}}) {
			my $cohort = 'FALSE';
			if(exists $data->{'pac'}{$name}{$pac}) {
				$cohort = 'TRUE';
			}
			$output .= "\t" . $cohort;
			
			my $is_hit = 'FALSE';
			if(exists $data->{'anno'}{$pac}{'2of3'}{$thres}{$name}) {
				if(1 == $data->{'anno'}{$pac}{'2of3'}{$thres}{$name}){
					$is_hit = 'TRUE';
				}
			}
			$output .= "\t" . $is_hit;

			my $is_hit2 = 'FALSE';
			if(exists $data->{'anno'}{$pac}{'hitlist'}{$name}) {
				if(1 == $data->{'anno'}{$pac}{'hitlist'}{$name}){
					$is_hit2 = 'TRUE';
				}
			}
			$output .= "\t" . $is_hit2;
			
			my $typeofhit = 'FALSE';
			if(exists $data->{'anno'}{$pac}{'hitlist'}{$name}) {
				$typeofhit = $data->{'anno'}{$pac}{'typeofhit'}{$name};
				
			}
			$output .= "\t" . $typeofhit;
			
			my $ratio = 0.000;
			if(exists $data->{'anno'}{$pac}{'ratio_patients_with_either'}{$name}) {
				$ratio = $data->{'anno'}{$pac}{'ratio_patients_with_either'}{$name};
			}
			$output .= "\t" . sprintf("%.4f", $ratio);
			
			my $vis = 0.000;
			if($is_hit eq 'TRUE') {
				$vis = $ratio;
			}
			$output .= "\t" . sprintf("%.4f", $vis);
			
			$ratio = 0.000;
			if(exists $data->{'anno'}{$pac}{'ratio_patients_with_snv'}{$name}) {
				$ratio = $data->{'anno'}{$pac}{'ratio_patients_with_snv'}{$name};
			}
			$output .= "\t" . sprintf("%.4f", $ratio);
			
			$vis = 0.000;
			if($is_hit eq 'TRUE') {
				$vis = $ratio;
			}
			$output .= "\t" . sprintf("%.4f", $vis);
			
			$ratio = 0.000;
			if(exists $data->{'anno'}{$pac}{'ratio_patients_with_cnv'}{$name}) {
				$ratio = $data->{'anno'}{$pac}{'ratio_patients_with_cnv'}{$name};
			}
			$output .= "\t" . sprintf("%.4f", $ratio);
			
			$vis = 0.000;
			if($is_hit eq 'TRUE') {
				$vis = $ratio;
			}
			$output .= "\t" . sprintf("%.4f", $vis);
		}
		$output .= "\n";
	}
	return $output;
}