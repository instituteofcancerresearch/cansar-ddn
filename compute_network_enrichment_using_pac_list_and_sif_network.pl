#!/usr/bin/perl
use warnings;
use strict;
use File::Temp qw/ tempfile /;
use Getopt::Long;


sub parse_hits_file($$);
sub parse_sifs_file($$);
sub compute_feature_nodes_enrichment_in_hits($);
sub generate_rscript();
sub generate_rscript_input($);
sub run_rscript($$);
sub import_R_stats($$);
sub generate_output($);


my $pac_list_file;
my $sif_file;
my $output_file;

GetOptions (
		#"help"			=> \$help,			#TO DO	
		"list=s"	=> \$pac_list_file,
		"sif=s" 	=> \$sif_file,
		"out=s"		=> \$output_file,
		);


die "Please provide a list of PAC identifiers to compute enrichment for\n" 
	unless($pac_list_file);
die "Please provide an interactome in sif format to compute enrichment with\n" 
	unless($sif_file);

my $data;

# PHASE 1: Basic data import
$data = parse_sifs_file($data, $sif_file);
$data = parse_hits_file($data, $pac_list_file);

# PHASE 3: Compute association enrichment
$data = compute_feature_nodes_enrichment_in_hits($data);

# PHASE 4: perform R statistics
my $rscript = generate_rscript();
my $rscript_input = generate_rscript_input($data);
my $rscript_output = run_rscript($rscript, $rscript_input);
$data = import_R_stats($data, $rscript_output);

# PHASE 5: generate output
my $output = generate_output($data);
if($output_file) {
	open(OUT, ">", $output_file) or die;
	print OUT $output;
	close(OUT);
}
else {
	print $output;
}

exit;



##########################

sub generate_output($){
	my $data = shift @_;
	my $output = 	
		"PAC\t" .
		"is_hit\t" .
		"Odds ratio\t" .
		"p-value\t" .
		"q-value\t" .
		"neighbours\t" .
		'hit/non hit ratio' . "\t" .
		"hits_and_1stN\t" .
		"hits_not_1stN\t" .
		"not_hits_1stN\t" .
		"not_hits_not_1stN\n";	
	foreach my $pac (keys  %{$data->{'enrichment'}}) {
		my $is_hit = 0;
		if(exists $data->{'hits'}{$pac}) {
			$is_hit = 1;
		}
		$output .=
		  	$pac . "\t" .
			$is_hit . "\t" .
			$data->{'enrichment'}{$pac}{'odds_r'} . "\t" .
			$data->{'enrichment'}{$pac}{'p-value'} . "\t" .
			$data->{'enrichment'}{$pac}{'q-value'} . "\t" .
			scalar(keys %{$data->{'enrichment'}{$pac}{'neighbour'}}) . "\t" .
			$data->{'enrichment'}{$pac}{'enrichment'} . "\t" .
			scalar(keys %{$data->{'enrichment'}{$pac}{'hits'}}) . "\t" .
			$data->{'enrichment'}{$pac}{'hits_not_neighbours'} . "\t" .
			$data->{'enrichment'}{$pac}{'not_hits_neighbours'} . "\t" .
			$data->{'enrichment'}{$pac}{'not_hits_not_neighbours'} . "\n";
	}
	return $output;
}

sub import_R_stats($$) {
	my ($data, $rscript_output_file) = @_;
	open(IN, $rscript_output_file) or die;
	my $header;
	while(defined(my $line = <IN>)) {
		unless($header) {
			$header++;
			next;
		}
		chomp($line);
		my ($pac, $pvalue, $odds_r, $odds_r_ll, $qvalue) = split(/\t/, $line);
		$data->{'enrichment'}{$pac}{'p-value'} = $pvalue;
		$data->{'enrichment'}{$pac}{'q-value'} = $qvalue;
		$data->{'enrichment'}{$pac}{'odds_r'} = $odds_r;
	}
	close(IN);
	return $data;
}

sub run_rscript($$) {
	my ($rscript, $input) = @_;
	my ($fh, $filename) = tempfile();
	my $command = "Rscript $rscript $input $filename";
	`$command`;
	return $filename;
}

sub generate_rscript_input($) {
	my $data = shift @_;
	my ($fh, $filename) = tempfile();
	foreach my $pac (keys  %{$data->{'enrichment'}})
	{
		print  $fh 	
			$pac . "\t" .
			scalar(keys %{$data->{'enrichment'}{$pac}{'hits'}}) . "\t" .
			$data->{'enrichment'}{$pac}{'hits_not_neighbours'} . "\t" .
			$data->{'enrichment'}{$pac}{'not_hits_neighbours'} . "\t" .
			$data->{'enrichment'}{$pac}{'not_hits_not_neighbours'} . "\n";
	}
	return $filename;
}


sub parse_hits_file($$) {
	my ($data, $file) = @_;
	open(IN, $file) or die "Could not read $file\n";
	while(defined(my $line = <IN>)) {
		chomp($line);		
		$data->{'hits'}{$line}++;
	}
	close(IN);
	return $data;
}

sub parse_sifs_file($$) {
	my ($data, $file) = @_;
	open(IN, $file) or die "Could not read $file\n";
	while(defined(my $line = <IN>))	{
		chomp($line);
		my ($left, $interaction, $right) = split(/\t/, $line);		
		# get rid of self interactions;
		next if($left eq $right);		
		$data->{'undirected_network'}{$left}{$right}++;
		$data->{'undirected_network'}{$right}{$left}++;
		#$data->{'directional_network'}{$left}{$right}{$interaction}++;
	}
	close(IN);
	return $data;
}


sub compute_feature_nodes_enrichment_in_hits($) {
	my $data = shift @_;	
	my $network_size = scalar(keys %{$data->{'undirected_network'}});
	my $total_hits = scalar(keys %{$data->{'hits'}});	
	foreach my $node (keys %{$data->{'undirected_network'}}) {
		my $hit_nodes = 0;
		my $first_neighbours = 0;
		my $is_hit_query = 0;
		if(exists $data->{'hits'}{$node}) {
			$is_hit_query = 1;
		}
		
		foreach my $interactor (keys %{$data->{'undirected_network'}{$node}}) {
			$first_neighbours++;
			$data->{'enrichment'}{$node}{'neighbour'}{$interactor}++;
			
			if(exists $data->{'hits'}{$interactor}) {
				$hit_nodes++;
				$data->{'enrichment'}{$node}{'hits'}{$interactor}++;
			}
		}
		my $enrichment = $hit_nodes / $first_neighbours;

		$data->{'enrichment'}{$node}{'enrichment'} = $enrichment;
		
		my $not_hits_neighbours = $first_neighbours -  $hit_nodes;
		my $hits_not_neighbours = $total_hits - $hit_nodes - $is_hit_query;
		my $not_hits_not_neighbours = $network_size - $total_hits - $not_hits_neighbours;
		$data->{'enrichment'}{$node}{'hits_not_neighbours'} = $hits_not_neighbours;
		$data->{'enrichment'}{$node}{'not_hits_neighbours'} = $not_hits_neighbours;
		$data->{'enrichment'}{$node}{'not_hits_not_neighbours'} = $not_hits_not_neighbours;
	}
	return $data;
}

sub generate_rscript() {
	my ($fh, $filename) = tempfile();
	print $fh 
		'#!/usr/bin/env Rscript' . "\n" .
		'require(qvalue)' . "\n" .
		'args = commandArgs(trailingOnly=TRUE)' . "\n" .
		'input <- read.table(as.character(args[1]), header=F, sep="\t", stringsAsFactors = FALSE)' . "\n" .
		'out <- as.character(args[2])' . "\n" .
		'fisher <- input[,c(2, 3, 4, 5)]' . "\n" .
		'p <- data.frame(t(apply(fisher, 1, function(x) {' . "\n" .
			'm <- matrix(as.numeric(x), ncol=2, byrow=T)' . "\n" .
    		'f <- fisher.test(m, alternative="greater")' . "\n" .
    			'return(c(#x,' . "\n" .
                		'p_val = f$p.value,' . "\n" .
                		'or = f$estimate[[1]],' . "\n" .
                		'or_ll = f$conf.int[[1]]))' . "\n" .
		'})))' . "\n" .
		'q <- qvalue(p$p_val)' . "\n" .
		'p <- cbind(input[,1], p, q_value = q$qvalues)' . "\n" .
		'colnames(p) <- c("PAC", "p-value", "OR", "OR_ll", "q-value")' . "\n" .
		'write.table(p, file=out, col.names = T, row.names = F, quote = F, sep = "\t")' . "\n";	
	return $filename;
}
