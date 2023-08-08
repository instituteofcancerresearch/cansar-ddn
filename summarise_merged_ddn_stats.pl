#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

sub read_best_druggability($$);
sub read_gene2pac($$);
sub read_annotation($$);
sub read_sif($$);
sub compute_drug_node_connectivity($);
sub print_drug_summary($);
sub print_hit_summary($);

my $gene2pac_file;
my $sif_file;
my $annotation_file;
my $best_druggability_file;
my $thres;
my $out_file;

GetOptions (
		#"help"				=> \$help,			#TO DO	
		"g2p=s"			=> \$gene2pac_file,
		"anno=s"		=> \$annotation_file,
		"drug=s" 		=> \$best_druggability_file,
		"sif=s" 		=> \$sif_file,
		"out=s"			=> \$out_file,
		);

my $data;
$data = read_gene2pac($data, $gene2pac_file);
$data = read_annotation($data, $annotation_file);
$data = read_best_druggability($data, $best_druggability_file);
$data = read_sif($data, $sif_file);
$data = compute_drug_node_connectivity($data);
my $drug_summary = print_drug_summary($data);
open(OUT, ">", $out_file . '.druggability.txt') or die;
print OUT $drug_summary;
close(OUT);
my $hit_summary = print_hit_summary($data);
open(OUT, ">", $out_file . '.hit_focus.txt') or die;
print OUT $hit_summary;
close(OUT);


sub read_gene2pac($$) {
	my ($data, $gene2pac_file) = @_;
	open(IN, $gene2pac_file) or die
		"Could not read $gene2pac_file\n";
	my $header;
	while(defined(my $line = <IN>))	{
		unless($header) {
			$header++;
			next;
		}
		chomp($line);
		my @cols = split(/\t/, $line);
		$data->{'p2g_map'}{$cols[1]} = $cols[0];
	}
	close(IN);
	return $data;	
}

sub read_annotation($$) {
	my ($data, $annotation_file) = @_;
	open(IN, $annotation_file) or die
		"Could not read $annotation_file\n";
	my $header;
	my %cohorts = ();
	while(defined(my $line = <IN>))	{
		chomp($line);
		my @cols = split(/\t/, $line);
		unless($header) {
			for(my $i  = 2; $i < scalar(@cols); $i = $i + 10) {
				my $cohort = $cols[$i];
				
				$cohort =~ s/cohort //;
				$cohorts{$cohort}{'node'} = $i;
				$cohorts{$cohort}{'hit'} = $i + 2;
				print STDERR "$cols[$i] \n";
			}
			$header++;
			next;
		}
		my $pac = $cols[0];
		my $gene = $cols[1];
		foreach my $cohort (keys %cohorts) {
			my $cohort_index = $cohorts{$cohort}{'node'};
			if('TRUE' eq $cols[$cohort_index]) {
				$data->{'cohorts'}{$cohort}{'node'}{$pac}{'gene'} =
					$data->{'p2g_map'}{$pac};
			}
			my $hit_index = $cohorts{$cohort}{'hit'};
			if('TRUE' eq $cols[$hit_index]) {
				$data->{'cohorts'}{$cohort}{'hit'}{$pac}++;
			}
		}
	}
	close(IN);
	return $data;
}

sub read_best_druggability($$) {
	my ($data, $best_druggability_file) = @_;
	open(IN, $best_druggability_file) or die
		"Could not read $best_druggability_file\n";
	my $header;
	while(defined(my $line = <IN>))	{
		unless($header) {
			$header++;
			next;
		}
		chomp($line);
		my ($pac, $status) = split(/\t/, $line);
		$data->{'best_druggability'}{$pac} = $status;
		foreach my $cohort (keys %{$data->{'cohorts'}}) {
			if(exists $data->{'cohorts'}{$cohort}{'node'}{$pac}) {
				$data->{'cohorts'}{$cohort}{'node'}{$pac}{'druggability'} =
					$status;
			}
		}
	}
	close(IN);	 
	return $data;	
}

sub read_sif($$) {
	my ($data, $sif_file) = @_;
	open(IN, $sif_file) or die
		"Could not read $sif_file\n";
	while(defined(my $line = <IN>))	 {
		chomp($line);
		my ($a, $type, $b) = split(/\t/, $line);
		foreach my $cohort (keys %{$data->{'cohorts'}}){
			if($data->{'cohorts'}{$cohort}{'node'}{$a}) {
				if($data->{'cohorts'}{$cohort}{'node'}{$b}) {
					$data->{'cohorts'}{$cohort}{'firstN'}{$a}{$b}++;
					$data->{'cohorts'}{$cohort}{'firstN'}{$b}{$a}++;
				}
			}
		}
	}
	close(IN);
	return $data;	
}

sub compute_drug_node_connectivity($) {
	my $data = shift @_;
	foreach my $cohort (keys %{$data->{'cohorts'}}){
		foreach my $node (keys %{$data->{'cohorts'}{$cohort}{'firstN'}}){
			my $status = $data->{'cohorts'}{$cohort}{'node'}{$node}{'druggability'};
			next unless($status);
			$data->{'cohorts'}{$cohort}{'drug'}{$node}{'status'} = $status;
			my $is_hit = 'FALSE';
			if(exists $data->{'cohorts'}{$cohort}{'hit'}{$node}) {
				$is_hit = 'TRUE';
			}
			foreach my $firstN (keys  %{$data->{'cohorts'}{$cohort}{'firstN'}{$node}}){
				if(exists $data->{'cohorts'}{$cohort}{'hit'}{$firstN}) {
					$data->{'cohorts'}{$cohort}{'drug_connect'}{$node}{'hits'}{$firstN}++;
				}
			}
		}
		
		foreach my $node (keys %{$data->{'cohorts'}{$cohort}{'hit'}}){
			my $status = $data->{'cohorts'}{$cohort}{'node'}{$node}{'druggability'};
			$status = 'challenging' unless($status);
			foreach my $firstN (keys  %{$data->{'cohorts'}{$cohort}{'firstN'}{$node}}){
				my $firstN_status = $data->{'cohorts'}{$cohort}{'node'}{$firstN}{'druggability'};
				next unless($firstN_status);
				$data->{'cohorts'}{$cohort}{'hit_connect'}{$node}{$firstN_status}{$firstN}++;
			}
		}
	}
	return $data;
}	

sub print_drug_summary($) {
	my $data = shift @_;
	my $output = "Gene\tPAC\tBest druggability";
	my @cohorts = sort keys %{$data->{'cohorts'}};
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' is_hit';
	}
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' number_of_1stN_hits';
	}
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' connected_hit_list';
	}
	$output .= "\n";
	my %pacs = ();
	foreach my $cohort (@cohorts) {
		foreach my $node (keys %{$data->{'cohorts'}{$cohort}{'drug_connect'}}) {
			$pacs{$node}++;
		}
	}
	foreach my $pac (sort keys %pacs) {
		my $gene = $data->{'p2g_map'}{$pac};
		my $status = $data->{'best_druggability'}{$pac};
		next unless($status);
		
		$output .= 	$gene . "\t" .
					$pac . "\t" .
					$status;
		
		foreach my $cohort (@cohorts) {
			my $is_hit = 'F';
			if(exists $data->{'cohorts'}{$cohort}{'hit'}{$pac}) {
				$is_hit = 'T';
			}
			$output .= 	"\t" . $is_hit;
		}
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'drug_connect'}{$pac}{'hits'}};
			$output .= 	"\t" . scalar(@firstN_pacs);
		}
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'drug_connect'}{$pac}{'hits'}};
			my %firstN_genes = ();
			foreach my $firstN_pac (@firstN_pacs) {
				my $gene =  $data->{'p2g_map'}{$firstN_pac};
				$firstN_genes{$gene}++;
			}
			$output .= 	"\t" .  join(", ", (sort keys %firstN_genes));
		}	
		$output	 .= "\n";			
	}
	return $output;
}


sub print_hit_summary($) {
	my $data = shift @_;
	my $output = "Gene\tPAC\tBest druggability";
	my @cohorts = sort keys %{$data->{'cohorts'}};
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' is_hit';
	}
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' number_of_Approved_1stN';
	}
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' list_of_Approved_1stN';
	}
	
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' number_of_Investigational_1stN';
	}
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' list_of_Investigational_1stN';
	}
	
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' number_of_Quality_probe_1stN';
	}
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' list_of_Quality_probe_1stN';
	}
	
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' number_of_3D_ligandable_1stN';
	}
	foreach my $cohort (@cohorts) {
		$output .= 	"\t" . $cohort . ' list_of_3D_ligandable_1stN';
	}
	$output .= "\n";
	my %pacs = ();
	foreach my $cohort (@cohorts) {
		foreach my $node (keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}}) {
			$pacs{$node}++;
		}
	}
	foreach my $pac (sort keys %pacs) {
		my $gene = $data->{'p2g_map'}{$pac};
		my $status = $data->{'best_druggability'}{$pac};
		$status = 'challenging' unless($status);
		
		$output .= 	$gene . "\t" .
					$pac . "\t" .
					$status;
		
		foreach my $cohort (@cohorts) {
			my $is_hit = 'F';
			if(exists $data->{'cohorts'}{$cohort}{'hit'}{$pac}) {
				$is_hit = 'T';
			}
			$output .= 	"\t" . $is_hit;
		}
		
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}{$pac}{'Approved'}};
			$output .= 	"\t" . scalar(@firstN_pacs);
		}
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}{$pac}{'Approved'}};
			my %firstN_genes = ();
			foreach my $firstN_pac (@firstN_pacs) {
				my $gene =  $data->{'p2g_map'}{$firstN_pac};
				$firstN_genes{$gene}++;
			}
			$output .= 	"\t" .  join(", ", (sort keys %firstN_genes));
		}
		
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}{$pac}{'Investigational'}};
			$output .= 	"\t" . scalar(@firstN_pacs);
		}
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}{$pac}{'Investigational'}};
			my %firstN_genes = ();
			foreach my $firstN_pac (@firstN_pacs) {
				my $gene =  $data->{'p2g_map'}{$firstN_pac};
				$firstN_genes{$gene}++;
			}
			$output .= 	"\t" .  join(", ", (sort keys %firstN_genes));
		}
		
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}{$pac}{'Quality_probe'}};
			$output .= 	"\t" . scalar(@firstN_pacs);
		}
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}{$pac}{'Quality_probe'}};
			my %firstN_genes = ();
			foreach my $firstN_pac (@firstN_pacs) {
				my $gene =  $data->{'p2g_map'}{$firstN_pac};
				$firstN_genes{$gene}++;
			}
			$output .= 	"\t" .  join(", ", (sort keys %firstN_genes));
		}
		
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}{$pac}{'3D_ligandable'}};
			$output .= 	"\t" . scalar(@firstN_pacs);
		}
		foreach my $cohort (@cohorts) {	
			my @firstN_pacs = keys %{$data->{'cohorts'}{$cohort}{'hit_connect'}{$pac}{'3D_ligandable'}};
			my %firstN_genes = ();
			foreach my $firstN_pac (@firstN_pacs) {
				my $gene =  $data->{'p2g_map'}{$firstN_pac};
				$firstN_genes{$gene}++;
			}
			$output .= 	"\t" .  join(", ", (sort keys %firstN_genes));
		}	
		$output	 .= "\n";			
	}
	return $output;
}	