#!/usr/bin/perl
use warnings;
no warnings 'recursion';
use strict;
use Getopt::Long;
use Graph::Undirected;
use Parallel::ForkManager;
use File::Temp qw/ tempfile tempdir /;
use diagnostics;

sub parse_list($);
sub add_special_nodes_to_hits($$);
sub generate_non_redundant_interactions($);
sub add_interactions_for_missing_hits($$);
sub add_special_nodes($$$);
sub add_trx_special_nodes($$$);
sub compute_enrichment($$$);
sub generate_interactome_object($$);
sub remove_non_hit_subnets($$);
sub compute_shortest_paths($$);
sub evaluate_shortest_paths($$$);
sub remove_non_hit_non_articulation($$$);
sub parse_druggability_file($$);
sub add_druggablility_to_graph($$$);
sub compute_bridges_and_simplify($$$);
sub remove_qp_3d_not_needed($$$);
sub generate_sif($$);
#sub find_bridges($$$);
sub parallel_compute_shortest_paths($$$);
sub compute_shortest_paths_for_a_hit($$$);
sub compute_shortest_paths_for_a_hit($$$);



my $sif_file; # = $ARGV[0];	# Non redundant interaction PAC_A $type PAC_B
my $hit_file; # = $ARGV[1];	# PACs that have been defined as hits
my $best_druggability_file;
my $special_nodes;
my $trx_special_nodes;
my $threads = 1;
my $out_file;
my $log_file;

GetOptions(
	"sif=s"		=>	\$sif_file, 
	"hits=s"	=>	\$hit_file,
	"drug=s"	=>	\$best_druggability_file,
	"special=s%"=>	\$special_nodes,
	"trx_special=s@" => \$trx_special_nodes,
	"threads=s" =>	\$threads,
	"out=s"		=>	\$out_file,
	"log=s"		=>	\$log_file,
);

$threads = 8 if($threads > 8);

die "SIF Hits and druggability files required\n"
	unless($sif_file && $hit_file && $best_druggability_file);

print STDERR "\n";
print STDERR localtime() . "\tSTAGE 1:\tParsing hits list file $hit_file\n";

my $hits = parse_list($hit_file);

if($special_nodes) {
	#$hits = add_special_nodes_to_hits($hits, $special_nodes);
}

print STDERR localtime() . "\tSTAGE 2:\tGenerating non redundant interactions\n";
my $nr_interactions = generate_non_redundant_interactions($sif_file);

# Need to add step for missing nodes
print STDERR localtime() . "\tSTAGE 3:\tAdding missing hits\n";
my $trimmed_nr_interactions = add_interactions_for_missing_hits($nr_interactions, $hits);

if($special_nodes) {
	print STDERR localtime() . "\tSTAGE 4:\tAdding special node and interactions\n";
	$trimmed_nr_interactions = add_special_nodes($nr_interactions, $trimmed_nr_interactions, $special_nodes);
}

if($trx_special_nodes) {
	print STDERR localtime() . "\tSTAGE 4:\tAdding special nodes and transcriptional interactions\n";
	$trimmed_nr_interactions = add_trx_special_nodes($nr_interactions, $trimmed_nr_interactions, $trx_special_nodes);
}

print STDERR localtime() . "\tSTAGE 5:\tComputing hit enrichments\n";
my $enrichments = compute_enrichment( $trimmed_nr_interactions, $hits, undef);

#print STDERR localtime() . "\tFinding viable non hit druggable candidates\n";
#my $eligible_druggable = parse_druggability_file($best_druggability_file, $enrichments);

print STDERR localtime() . "\tSTAGE 6:\tGenerating interactome and annotating edge type\n";
my $nr_interactome = generate_interactome_object($trimmed_nr_interactions, undef);
$nr_interactome->biconnectivity_clear_cache;
if($nr_interactome->is_connected) {
	print STDERR localtime() . "\tSTAGE 6a:\tConnected interactome after generate_interactome_object\n";
}
else {
	print STDERR localtime() . "\tSTAGE 6a:\tTrying to remove subgraphs containing non hits\n";
	$nr_interactome = remove_non_hit_subnets($nr_interactome, $hits);
}

print STDERR localtime() . "\tSTAGE 7:\tComputing shortest paths between mapped hits\n";
#my $shortest = compute_shortest_paths($nr_interactome, $hits);
my $shortest = parallel_compute_shortest_paths($nr_interactome, $hits, $threads);


print STDERR localtime() . "\tSTAGE 8:\tScoring shortest paths for length and enrichment\n";
my $best = evaluate_shortest_paths($shortest, $hits, $enrichments);

# sort out articulation points
print STDERR localtime() . "\tSTAGE 9:\tRemoving non-hit nodes unless they are articulation points\n";
my $trimmed_subnet = remove_non_hit_non_articulation($best, $nr_interactome, $hits);
$trimmed_subnet->biconnectivity_clear_cache;
if($trimmed_subnet->is_connected) {
	print STDERR localtime() . "\tSTAGE 9:\tConnected graph after remove_non_hit_non_articulation\n";
}
else {
	my @cc = $trimmed_subnet->connected_components();
	print STDERR localtime() . "\tSTAGE 9:\tgraph components: " . scalar( @cc ) . "\n";
	print STDERR localtime() . "\tSTAGE 9:\tTrying to remove subgraphs containing non hits\n";
	$trimmed_subnet = remove_non_hit_subnets($trimmed_subnet, $hits);
}


print STDERR localtime() . "\tSTAGE 10:\tFinding viable non hit druggable candidates\n";
my $eligible_druggable = parse_druggability_file($best_druggability_file, $enrichments);



print STDERR localtime() . "\tSTAGE 11:\tAdding eligible druggability hits to graph and subsetting\n";
my $druggable_subnet = add_druggablility_to_graph($trimmed_subnet, $nr_interactome, $eligible_druggable);
$druggable_subnet->biconnectivity_clear_cache;
if($druggable_subnet->is_connected) {
	print STDERR localtime() . "\tConnected graph after add_druggablility_to_graph\n";
}
else {
	my @cc = $druggable_subnet->connected_components();
	print STDERR localtime() . "\tgraph components: " . scalar( @cc ) . "\n";
	$druggable_subnet = remove_non_hit_subnets($druggable_subnet, $hits);
}


print STDERR localtime() . "\tSTAGE 12:\tSimplifying network\n";
my $trimmed_druggable_subnet = compute_bridges_and_simplify($druggable_subnet, $eligible_druggable, $hits);
$trimmed_druggable_subnet->biconnectivity_clear_cache;
if($trimmed_druggable_subnet->is_connected) {
	print STDERR localtime() . "\tConnected graph after compute_bridges_and_simplify\n";
}
else {
	my @cc = $trimmed_druggable_subnet->connected_components();
	print STDERR localtime() . "\tgraph components: " . scalar( @cc ) . "\n";
	$trimmed_druggable_subnet = remove_non_hit_subnets($trimmed_druggable_subnet, $hits);

}



print STDERR localtime() . "\tSTAGE 13:\tFurther Simplifying network\n";
my $further_trimmed_druggable = remove_qp_3d_not_needed($trimmed_druggable_subnet,$eligible_druggable, $hits);
$further_trimmed_druggable->biconnectivity_clear_cache;
if($further_trimmed_druggable->is_connected) {
	print STDERR localtime() . "\tConnected graph after remove_qp_3d_not_needed\n";
}
else {
	my @cc = $further_trimmed_druggable->connected_components();
	print STDERR localtime() . "\tgraph components: " . scalar( @cc ) . "\n";
	$further_trimmed_druggable = remove_non_hit_subnets($further_trimmed_druggable, $hits);
}


#if($further_trimmed_druggable) {
#	print "The graph is $further_trimmed_druggable\n";
#}


my $sif_output = generate_sif($further_trimmed_druggable, $trimmed_nr_interactions);

my $fail;
if($out_file) {
	open(OUT, ">", $out_file) or $fail++;
	unless($fail) {
		print STDERR localtime() . "\tPrinting SIF network to file $out_file\n";
		print OUT $sif_output;
		close(OUT);
	}
}
else {
	print $sif_output;
} 
if($fail) {
	print $sif_output;
}							

print STDERR localtime() . "\tFinished\n";

exit;


###
###
###	SUBROUTINES
###
###

sub parse_list($)
{
	my $file = shift @_;
	open(IN, $file) or die;
	my %data = ();
	while(defined(my $line = <IN>))
	{
		chomp($line);
		$data{$line}++;
	}
	close(IN);
	print STDERR localtime() . "\tParsed " . scalar(keys %data) . " hits\n";
	return \%data;
}

#
#	add special nodes to hit list if defined
#
sub add_special_nodes_to_hits($$) {
	my ($hits, $special_nodes) = @_;
	my $old = scalar(keys %$hits);
	my $special_node = $special_nodes->{'node'};
	print STDERR localtime() . "\tAdding special_node $special_node to hit list\n";
	$hits->{$special_node}++;
	my $new = scalar(keys %$hits);
	print STDERR localtime() . "\t$old --> $new hits\n";
	return $hits;
}

#
#	Read in interactions into a hashmap and makes them non redundant
#
sub generate_non_redundant_interactions($) {
	my $file = shift @_;	
	my %input = ();
	my %seen = ();
	my %output = ();
	
	# first parse interactome file
	# This is a simple sif file containing the right type interactions
	# Interactors have PAC accessions
	open(IN, $file) or die
		"Could not access file $file\n";
	while(defined(my $line = <IN>)) {
		chomp($line);
		my ($left, $type, $right) = split(/\t/, $line);
		next if($left eq $right);
		$input{$type}{$left}{$right}++;
	}
	close(IN);
	
	# generate non redundant	
	foreach my $left (keys %{$input{'OncoKB->'}}) {
		foreach my $right  (keys %{$input{'OncoKB->'}{$left}}) {
			$output{'OncoKB->'}{$left}{$right}++;
			$seen{$left}{$right}++;
		}
	}
	foreach my $left (keys %{$input{'Signor->'}}) {
		foreach my $right  (keys %{$input{'Signor->'}{$left}}) {
			next if(exists $seen{$left}{$right});
			$output{'Signor->'}{$left}{$right}++;
			$seen{$left}{$right}++;
		}
	}
	foreach my $left (keys %{$input{'Signor-|'}}) {
		foreach my $right  (keys %{$input{'Signor-|'}{$left}}) {
			next if(exists $seen{$left}{$right});
			$output{'Signor-|'}{$left}{$right}++;
			$seen{$left}{$right}++;
		}
	}
	foreach my $left (keys %{$input{'reaction->'}}) {
		foreach my $right  (keys %{$input{'reaction->'}{$left}}) {
			next if(exists $seen{$left}{$right});
			$output{'reaction->'}{$left}{$right}++;
			$seen{$left}{$right}++;
		}
	}
	foreach my $left (keys %{$input{'direct (LigInt)'}}) {
		foreach my $right  (keys %{$input{'direct (LigInt)'}{$left}}) {
			next if(exists $seen{$left}{$right});
			next if(exists $seen{$right}{$left});
			$output{'direct (LigInt)'}{$left}{$right}++;
			$seen{$left}{$right}++;
			$seen{$right}{$left}++;
		}
	}
	foreach my $left (keys %{$input{'direct (X-ray)'}}) {
		foreach my $right  (keys %{$input{'direct (X-ray)'}{$left}}) {
			next if(exists $seen{$left}{$right});
			next if(exists $seen{$right}{$left});
			$output{'direct (X-ray)'}{$left}{$right}++;
			$seen{$left}{$right}++;
			$seen{$right}{$left}++;
		}
	}
	foreach my $left (keys %{$input{'direct (nonPDB)'}}) {
		foreach my $right  (keys %{$input{'direct (nonPDB)'}{$left}}) {
			next if(exists $seen{$left}{$right});
			next if(exists $seen{$right}{$left});
			$output{'direct (nonPDB)'}{$left}{$right}++;
			$seen{$left}{$right}++;
			$seen{$right}{$left}++;
		}
	}
	foreach my $left (keys %{$input{'transcriptional'}}) {
		foreach my $right  (keys %{$input{'transcriptional'}{$left}}) {
			next if(exists $seen{$left}{$right});
			next if(exists $seen{$right}{$left});
			$output{'transcriptional'}{$left}{$right}++;
			$seen{$left}{$right}++;
		}
	}
	foreach my $left (keys %{$input{'complex'}}) {
		foreach my $right  (keys %{$input{'complex'}{$left}}) {
			next if(exists $seen{$left}{$right});
			next if(exists $seen{$right}{$left});
			$output{'complex'}{$left}{$right}++;
			$seen{$left}{$right}++;
			$seen{$right}{$left}++;
		}
	}	
	my $total_interactor_counter = 0;
	foreach my $type (sort keys %output) {
		my $interaction_counter = 0;
		foreach my $nodeA (sort keys %{$output{$type}}) {
			foreach my $nodeB  (sort keys %{$output{$type}{$nodeA}}) {
				$interaction_counter++;
			}
		}
		$total_interactor_counter += $interaction_counter;
		print STDERR localtime() . "\tAdded $interaction_counter $type interactions\n";
	}	
	print STDERR localtime() . "\tAdded $total_interactor_counter interactions\n";	
	return \%output;
}

#
#	Maps hits to canonical interactome and subsequently tries to add
#	missing hits via trx and complex interactions using enrichment
#
sub add_interactions_for_missing_hits($$) {
	my ($interactions, $hits) = @_;
	my %output = ();
	my %seen_hits = ();
	my %strict_nodes = ();
	# First add standard interaction types
	foreach my $type (keys %$interactions) {
		next if($type eq 'complex');
		next if($type eq 'complex_low');
		next if($type eq 'transcriptional');
		foreach my $left (keys %{$interactions->{$type}}) {
			foreach my $right (keys %{$interactions->{$type}{$left}}) {
				if(exists $hits->{$left}) {
					$seen_hits{$left}++;
				}
				if(exists $hits->{$right}) {
					$seen_hits{$right}++;
				}
				$output{$type}{$left}{$right}++;
				$strict_nodes{$left}++;
				$strict_nodes{$right}++;
			}
		}		
	}
	print STDERR  localtime() . "\t" . scalar(keys %seen_hits) . " hits in strict interactome\n";
	# Find which are the missing hits	
	my %missing_hits = ();
	foreach my $hit (keys %$hits) {
		next if(exists $seen_hits{$hit});
		$missing_hits{$hit}++;
	}
	print STDERR localtime() . "\tWill try to add " . scalar(keys %missing_hits) . " missing hits\n";	
	
	# Add all direct trx associations
	my %added_hits = ();
	foreach my $hit (keys %$hits) {
		next if(exists $seen_hits{$hit});
		foreach my $seen_hit (keys %seen_hits) {
			if(exists $interactions->{'transcriptional'}{$hit}{$seen_hit}) {
				$output{'transcriptional'}{$hit}{$seen_hit}++;
				$seen_hits{$hit}++;
				$added_hits{'direct'}{$hit} = 'transcriptional';
			}
			elsif(exists $interactions->{'transcriptional'}{$seen_hit}{$hit}) {
				$output{'transcriptional'}{$seen_hit}{$hit}++;
				$seen_hits{$hit}++;
				$added_hits{'direct'}{$hit} = 'transcriptional';
			}
		}
	}
	
	# Add all direct complex associations
	foreach my $hit (keys %$hits) {
		next if(exists $seen_hits{$hit});
		foreach my $seen_hit (keys %seen_hits) {
			if(exists $interactions->{'complex'}{$hit}{$seen_hit}) {
				$output{'complex'}{$hit}{$seen_hit}++;
				$seen_hits{$hit}++;
				$added_hits{'direct'}{$hit} = 'complex';
			}
			elsif(exists $interactions->{'complex'}{$seen_hit}{$hit}) {
				$output{'complex'}{$seen_hit}{$hit}++;
				$seen_hits{$hit}++;
				$added_hits{'direct'}{$hit} = 'complex';
			}
		}
	}
	
	# Now try to add remaining missing hits by casting the net wider:
	my $preliminary_enrichment = compute_enrichment(\%output, $hits, undef);
	
	foreach my $hit (keys %$hits) {
		next if(exists $seen_hits{$hit});
		my $best_enrichment;
		my $best_partner;
		
		foreach my $right (keys %{$interactions->{'transcriptional'}{$hit}}) {
			next unless(exists $strict_nodes{$right});
			my $enrichment = $preliminary_enrichment->{$right}{'ratio'};
			if($best_partner) {
				if($enrichment > $best_enrichment) {
					$best_partner = $right;
					$best_enrichment = $enrichment;
				}
			}
			else {
				$best_partner = $right;
				$best_enrichment = $enrichment;
			}
		}
		foreach my $left (keys %{$interactions->{'transcriptional'}}) {
			next unless(exists $strict_nodes{$left});
			next unless(exists $interactions->{'transcriptional'}{$left}{$hit});
			my $enrichment = $preliminary_enrichment->{$left}{'ratio'};
			if($best_partner) {
				if($enrichment > $best_enrichment) {
					$best_partner = $left;
					$best_enrichment = $enrichment;
				}
			}
			else {
				$best_partner = $left;
				$best_enrichment = $enrichment;
			}
		}
		if($best_partner) {
			$output{'transcriptional'}{$hit}{$best_partner}++;
			$seen_hits{$hit}++;
			$added_hits{'indirect'}{$hit} = 'transcriptional';
		}	
	}
	# Recompute enrichment and try adding complex interactions:
	$preliminary_enrichment = compute_enrichment(\%output, $hits, undef);
	foreach my $hit (keys %$hits) {
		next if(exists $seen_hits{$hit});
		my $best_enrichment;
		my $best_partner;
		foreach my $right (keys %{$interactions->{'complex'}{$hit}}) {
			next unless(exists $strict_nodes{$right});
			my $enrichment = $preliminary_enrichment->{$right}{'ratio'};
			if($best_partner) {
				if($enrichment > $best_enrichment) {
					$best_partner = $right;
					$best_enrichment = $enrichment;
				}
			}
			else {
				$best_partner = $right;
				$best_enrichment = $enrichment;
			}
		}
		foreach my $left (keys %{$interactions->{'complex'}}) {
			next unless(exists $strict_nodes{$left});
			next unless(exists $interactions->{'complex'}{$left}{$hit});
			my $enrichment = $preliminary_enrichment->{$left}{'ratio'};
			if($best_partner) {
				if($enrichment > $best_enrichment) {
					$best_partner = $left;
					$best_enrichment = $enrichment;
				}
			}
			else {
				$best_partner = $left;
				$best_enrichment = $enrichment;
			}
		}
		if($best_partner) {
			$output{'complex'}{$hit}{$best_partner}++;
			$seen_hits{$hit}++;
			$added_hits{'indirect'}{$hit} = 'complex';
		}	
	}
	%missing_hits = ();
	foreach my $hit (keys %$hits) {
		next if(exists $seen_hits{$hit});
		$missing_hits{$hit}++;
	}
	
	foreach my $type (sort keys %added_hits) {
		foreach my $hit (sort keys %{$added_hits{$type}}) {
			my $edge = $added_hits{$type}{$hit};
			print STDERR localtime() . "\tAdded $hit --> $type association with hits of class $edge\n";
		}
	}
		
	my %final_hits = ();
	foreach my $class (sort keys %output) {
		foreach my $left (keys %{$output{$class}}) {
			if(exists $hits->{$left}) {
				$final_hits{$left}++;
			}
			foreach my $right (keys %{$output{$class}{$left}}) {
				if(exists $hits->{$right}) {
					$final_hits{$right}++;
				}
			}
		}
	}
	print STDERR localtime() . "\t" . scalar(keys %final_hits) . " of " . scalar(keys %$hits) . " hits can be investigated\n";
	print STDERR localtime() . "\t" . scalar(keys %missing_hits) . 
			" hits cannot be used: " . join(", ", (sort keys %missing_hits)) . "\n";
	return \%output;	
}

#
#	Adds user specified interaction types for a user defined node
#	e.g. transcriptional interactions for ESR1
#
sub add_special_nodes($$$) {
	my ($nr_interactions, $trimmed_nr_interactions, $special_nodes) = @_;
	my $counter = 0;
	my $special_node = $special_nodes->{'node'};
	my $special_interaction = $special_nodes->{'type'};
	
	# check which interactions we already have.  
	# We are only going to add missing interactions
	my %existsing_interactions = ();
	foreach my $type (keys %$trimmed_nr_interactions) {
		foreach my $left (keys %{$trimmed_nr_interactions->{$type}}) {
			foreach my $right (keys %{$trimmed_nr_interactions->{$type}{$left}}) {
				$existsing_interactions{$left}{$right}++;
				$existsing_interactions{$right}{$left}++;
			}
		}	
	}
	# Now grep interactions to add:
	# They should contain the special node and
	# the connecting edge should not exist in
	# the trimmed network
	foreach my $left (keys %{$nr_interactions->{$special_interaction}}) {
		foreach my $right (keys %{$nr_interactions->{$special_interaction}{$left}}) {
			next if($left eq $right); # dont want self interactions
			if($left eq $special_node || $right eq $special_node) {
				unless(exists $existsing_interactions{$left}{$right}) {
					$trimmed_nr_interactions->{$special_interaction}{$left}{$right}++;
					$counter++;
				}	
			}
		}
	}
	print STDERR localtime() . "\tAdded $counter $special_interaction interactions for " .
			" special node $special_node\n";		
	return $trimmed_nr_interactions;
}

#
#	Adds user specified interaction types for a user defined node
#	e.g. transcriptional interactions for ESR1
#
sub add_trx_special_nodes($$$) {
	my ($nr_interactions, $trimmed_nr_interactions, $special_nodes) = @_;
	my $counter = 0;
	my $special_interaction = 'transcriptional';
	
	# check which interactions we already have.  
	# We are only going to add missing interactions
	my %existsing_interactions = ();
	foreach my $type (keys %$trimmed_nr_interactions) {
		foreach my $left (keys %{$trimmed_nr_interactions->{$type}}) {
			foreach my $right (keys %{$trimmed_nr_interactions->{$type}{$left}}) {
				$existsing_interactions{$left}{$right}++;
				$existsing_interactions{$right}{$left}++;
			}
		}	
	}
	# Now grep interactions to add:
	# They should contain the special node and
	# the connecting edge should not exist in
	# the trimmed network
	foreach my $left (keys %{$nr_interactions->{$special_interaction}}) {
		foreach my $right (keys %{$nr_interactions->{$special_interaction}{$left}}) {
			next if($left eq $right); # dont want self interactions
			foreach my $special_node (@$special_nodes) {
				if($left eq $special_node || $right eq $special_node) {
					unless(exists $existsing_interactions{$left}{$right}) {
						$trimmed_nr_interactions->{$special_interaction}{$left}{$right}++;
						$counter++;
					}	
				}	
			}
		}
	}
	print STDERR localtime() . "\tAdded $counter $special_interaction interactions for " .
			scalar(@$special_nodes) . " special nodes\n";		
	return $trimmed_nr_interactions;
}

#
#	Establish hit ratio for first neighbours of each interactome node
#
sub compute_enrichment($$$) {
	my ($interactions, $hits, $stingent) = @_;
	my %enrichment = ();
	foreach my $type (keys %$interactions) {
		if($stingent) {
			next if($type eq 'complex');
			next if($type eq 'complex_low');
			next if($type eq 'transcriptional');
		}
		foreach my $left (keys %{$interactions->{$type}}) {
			foreach my $right (keys %{$interactions->{$type}{$left}}) {
				if(exists $hits->{$right}) {
					$enrichment{$left}{'hits'}{$right}++;
				}
				$enrichment{$left}{'neighbours'}{$right}++;
				if(exists $hits->{$left}) {
					$enrichment{$right}{'hits'}{$left}++;
				}
				$enrichment{$right}{'neighbours'}{$left}++;
			}
		}
	}
	foreach my $node (keys %enrichment) {
		my $hits = scalar( keys %{$enrichment{$node}{'hits'}});
		my $neighbours = scalar( keys %{$enrichment{$node}{'neighbours'}});
		if($neighbours) {
			$enrichment{$node}{'ratio'} = $hits / $neighbours;
		}
		else {
			$enrichment{$node}{'ratio'} = 0;
		}
		unless($enrichment{$node}{'hits'}) {
			$enrichment{$node}{'hits'} = 0;
		}
	}
	return \%enrichment;	
}

#
#	Converts interactions hashmap to interactome object
#
sub generate_interactome_object($$) {
	my ($interactions, $stingent) = @_;
	my $g = Graph::Undirected->new;
	my $edge_counter = 0;
	foreach my $type (keys %$interactions) {
		if($stingent) {
			next if($type eq 'complex');
			next if($type eq 'complex_low');
			next if($type eq 'transcriptional');
		}
		foreach my $left (sort keys %{$interactions->{$type}}) {
			foreach my $right (sort keys %{$interactions->{$type}{$left}}) {
				$g->set_edge_attribute($left, $right, 'type', $type);
				$edge_counter++;
			}
		}
	}
	my @fnodes = $g->vertices;
	my $hit_counter = 0;
	foreach my $node (@fnodes) {
		if(exists $hits->{$node}) {
			$hit_counter++;
		}
	}	
	my @fedges = $g->edges;
	print STDERR localtime() . "\tGenerated network with " . scalar(@fnodes) . " nodes and " .
						scalar(@fedges) . " edges containing $hit_counter hits\n";		

					
	return $g;
}

sub remove_non_hit_subnets($$) {
	my ($graph, $hits) = @_;
	$graph->biconnectivity_clear_cache;
	my @nodes = $graph->vertices;
	my @edges = $graph->edges;
	print STDERR localtime() . "\tInput graph contains " . scalar(@nodes) . " nodes and " .
							scalar(@edges) . " edges\n";
	if($graph->is_connected) {
		return $graph;
	}
	else {
		my @ccs = $nr_interactome->connected_components();
		print STDERR localtime() . "\tNetwork contains " . scalar(@ccs) . " subnets\n";
		my $subnet_counter = 0;
		foreach my $cc (@ccs) {
			$subnet_counter++;
			my $subnet_size = scalar(@$cc);
			my %subnet_hits = ();
			foreach my $subnet_node (@$cc) {
				if(exists $hits->{$subnet_node}) {
					$subnet_hits{$subnet_node}++;
				}
			}
			if(scalar(keys %subnet_hits) > 0) {
				print STDERR localtime() . "\tRetaining subnet $subnet_counter of size $subnet_size as it contains " . 
					scalar(keys %subnet_hits) . " hits\n";
				print STDERR localtime() . "\thit nodes: " . join(" ", (sort keys %subnet_hits)) . "\n";
			}
			else {
				foreach my $subnet_node (@$cc) {
					#print STDERR localtime() . "\tDeleting node $subnet_node in subnet $subnet_counter\n";
					$graph->delete_vertex($subnet_node);
				}
				print STDERR localtime() . "\tDeleted subnet $subnet_counter of size $subnet_size as it contained no hits\n";
			}
		}
	}
	@nodes = $graph->vertices;	
	@edges = $graph->edges;
	print STDERR localtime() . "\tOutput graph contains " . scalar(@nodes) . " nodes and " .
							scalar(@edges) . " edges\n";
	return $graph;
}


sub parallel_compute_shortest_paths($$$) {
	my ($graph, $hits, $threads) = @_;
	my @graph_vertices = $graph->vertices;
	my %shortest = ();
	$graph->SPT_Dijkstra_clear_cache;
	my %hash = map { $_ => 1 } @graph_vertices;
	# Establish which hits are in the interactome:
	foreach my $hit (keys %$hits) {
		if(exists $hash{$hit}) {
			$shortest{'present'}{$hit}++;
		}
		else {
			$shortest{'missing'}{$hit}++;
		}
	}
	
	my %temp_files = ();
	# Now compute shortest paths for the present nodes;
	my @present_nodes = sort keys %{$shortest{'present'}};
	$graph->SPT_Dijkstra_clear_cache;
	my $counter_down = scalar(@present_nodes);
	print STDERR localtime() . "\tLeft: ";
	
	my $pm = new Parallel::ForkManager($threads);
	
	for(my $u = 0; $u < scalar(@present_nodes); $u++) {
		print STDERR $counter_down . " ";
		$counter_down--;
		my $temp_file = File::Temp->new();
		$temp_files{$temp_file}++;
		#print STDERR localtime() . "\tComputing in between hits shortest paths for " .
		#		$present_nodes[$u] . "\tleft: $counter_down\n";
		my $pid = $pm->start and next;	
		my $results = compute_shortest_paths_for_a_hit($graph, $hits, $present_nodes[$u]);
		open(OUT, ">", $temp_file) or die;
		# write to temp file
		foreach my $a (sort keys %{$results->{'shortest'}}) {
			foreach my $length (sort keys %{$results->{'shortest'}{$a}}) {
				foreach my $b (sort keys %{$results->{'shortest'}{$a}{$length}}) {
					my $path = $results->{'shortest'}{$a}{$length}{$b};
					#$shortest{'shortest'}{$a}{$length}{$b} = $path;
					print OUT $a . "\t" . $length . "\t" . $b . "\t" . join(",", @$path) . "\n";
					
				}
			}
		}
		close(OUT);
		$pm->finish;	
	}
	$pm->wait_all_children;
	print STDERR "\n";
	my $line_counter = 0;
	foreach my $temp_file (keys %temp_files) {
		open(IN, $temp_file) or die;
		while(defined(my $line = <IN>)) {
			chomp($line);
			my ($a, $length, $b, $path) = split(/\t/, $line);
			my @paths = split(",", $path);
			#####print "Unfiltered shortests: $length\t" . $a . " --> " . join(" --> ", @paths) . " --> " . $b . "\n";
			$shortest{'shortest'}{$a}{$length}{$b} = \@paths;
			$line_counter++;
		}
		close(IN);
	}
	my %hit_counter = ();
	foreach my $a (keys %{$shortest{'shortest'}}) {
		foreach my $length (keys %{$shortest{'shortest'}{$a}}) {
			foreach my $b (keys %{$shortest{'shortest'}{$a}{$length}}) {
				if(exists $hits->{$a}) {
					$hit_counter{$a}++;
					if(exists $hits->{$b}) {
						$shortest{'connectivity'}{$a}{$b}++;
						$hit_counter{$b}++;
					}	
				}
			}
		}
	}
	foreach my $a (sort keys %{$shortest{'connectivity'}}) {
		print STDERR  localtime() . "\t$a\t" . 
					scalar(keys %{$shortest{'connectivity'}{$a}} ) . "\t" . 
					scalar(keys %hit_counter) . "\n";
	}
	print STDERR localtime() . "\tShortest paths contain " . scalar(keys %hit_counter) . " hits\n";	
	
				
	return \%shortest;

}

sub compute_shortest_paths_for_a_hit($$$) {
	my ($graph, $hits, $node) = @_;
	my @graph_vertices = $graph->vertices;
	my %shortest = ();
	$graph->SPT_Dijkstra_clear_cache;
	my %hash = map { $_ => 1 } @graph_vertices;
	# Establish which hits are in the interactome:
	foreach my $hit (sort keys %$hits) {
		if(exists $hash{$hit}) {
			$shortest{'present'}{$hit}++;
		}
		else {
			$shortest{'missing'}{$hit}++;
		}
	}
	my @present_nodes = sort keys %{$shortest{'present'}};
	for(my $v = 0; $v < scalar(@present_nodes); $v++) {
		next if($node eq $present_nodes[$v]);
		my $results = compute_shortest_paths_for_single_node_pair($graph, $node, $present_nodes[$v]);
		foreach my $a (sort keys %$results) {
			foreach my $length (sort keys %{$results->{$a}}) {
				foreach my $b (sort keys %{$results->{$a}{$length}}) {
					my $path = $results->{$a}{$length}{$b};
					$shortest{'shortest'}{$a}{$length}{$b} = $path;
				}
			}
		}
		
	}
	return \%shortest;
}

sub compute_shortest_paths_for_single_node_pair($$$) {
	my ($graph, $nodeA, $nodeB) = @_;
	my @graph_vertices = $graph->vertices;
	my %shortest = ();	
	my @path = $graph->SP_Dijkstra($nodeA, $nodeB);
	my $length = scalar(@path);
	$shortest{$nodeA}{$length}{$nodeB} = \@path;
	return \%shortest;
}


#
#	Identifies all shortest paths between the mappable hits in the interactome
#	Replaced by the parallel thread processing routines
sub compute_shortest_paths($$) {
	my ($graph, $hits) = @_;
	my @graph_vertices = $graph->vertices;
	my %shortest = ();
	$graph->SPT_Dijkstra_clear_cache;
	my %hash = map { $_ => 1 } @graph_vertices;
	# Establish which hits are in the interactome:
	foreach my $hit (keys %$hits) {
		if(exists $hash{$hit}) {
			$shortest{'present'}{$hit}++;
		}
		else {
			$shortest{'missing'}{$hit}++;
		}
	}
	#print STDERR localtime() . "\tMissing nodes: " . join(" ", sort keys %{$shortest{'missing'}}) . "\n";
	#print STDERR localtime() . "\tComputing shortest paths between " . scalar(keys %{$shortest{'present'}}) . 
	#		" of " . scalar(keys %$hits) . " hits in the network\n";
			
	# Now compute shortest paths for the present nodes;
	my @present_nodes = sort keys %{$shortest{'present'}};
	$graph->SPT_Dijkstra_clear_cache;
	my $counter_down = scalar(@present_nodes);
	print STDERR localtime() . "\tLeft: ";
	for(my $u = 0; $u < scalar(@present_nodes); $u++) {
	print STDERR $counter_down . " ";
	$counter_down--;
		#print STDERR localtime() . "\tComputing in between hits shortest paths for " .
		#	$present_nodes[$u] . "\tleft: $counter_down\n";
			
		for( my $v = $u; $v < scalar(@present_nodes); $v++) {
			next if($u eq $v);
			
			my @path = $graph->SP_Dijkstra($present_nodes[$u], $present_nodes[$v]);
			my $length = scalar(@path);
			$shortest{'shortest'}{$present_nodes[$u]}{$length}{$present_nodes[$v]} = \@path;
		}
	}
	print STDERR "\n";
	return \%shortest;
}

#
#	Select best of shortest paths between the mappable hits in the interactome
#	Best: shortest path of highest enrichment
#
sub evaluate_shortest_paths($$$) {
	my ($shortest, $hits, $enrichments) = @_;
	my %visits = ();
	my %best = ();
	my %reconfigured = ();
	foreach my $hitA (sort keys %{$shortest->{'shortest'}}) {
		foreach my $length (keys %{$shortest->{'shortest'}{$hitA}}) {
			foreach my $hitB (sort keys %{$shortest->{'shortest'}{$hitA}{$length}}) {
				my $path = $shortest->{'shortest'}{$hitA}{$length}{$hitB};
				$reconfigured{$hitA}{$hitB}{$length} = $path;
			}
		}
	}
	
	# Compute how often a node is visited
	#print STDERR localtime() . "\t\tComputing node visit frequency using shortest paths\n";
	foreach my $hitA (sort keys %reconfigured) {
		foreach my $hitB (sort keys %{$reconfigured{$hitA}}) {
			my @lengths = sort {$a <=> $b} keys %{$reconfigured{$hitA}{$hitB}};
			my $min_length = $lengths[0];
			my $start_neigh = 1;
			my $end_neigh = $min_length - 1;	
			my $path = $reconfigured{$hitA}{$hitB}{$min_length};
			foreach my $node (@$path[$start_neigh..$end_neigh]) {
				$visits{$node}++;
			}
		}
	}
	

	# Now evaluate same length paths
	#print STDERR localtime() . "\t\tScoring paths\n";
	foreach my $hitA (sort keys %reconfigured) {
		foreach my $hitB (sort keys %{$reconfigured{$hitA}}) {
			my @lengths = sort {$a <=> $b} keys %{$reconfigured{$hitA}{$hitB}};
			my $min_length = $lengths[0];
			my $start_neigh = 1;
			my $end_neigh = $min_length - 1;	
			my $path = $reconfigured{$hitA}{$hitB}{$min_length};
			my $score = 0;
			my $enrichment = 1;
			foreach my $node (@$path[$start_neigh..$end_neigh]) {
				$score += 50 if(exists $hits->{$node});	
				if(exists $visits{$node}) {
					$score += $visits{$node};
				}
				if(exists $enrichments->{$node}{'enrichment'}) {
					$enrichment *= $enrichments->{$node}{'enrichment'};
				}
				else {
					unless(exists $hits->{$node}) {
						$enrichment *= 0.001; 
					}
				}
			}
			if(exists $best{$hitA}{'enrichment'}) {
				my $prior_enrichment = $best{$hitA}{'enrichment'};
				if($enrichment > $prior_enrichment) {
					$best{$hitA}{$hitB}{'enrichment'} = $enrichment;
					$best{$hitA}{$hitB}{'score'} = $score;
					$best{$hitA}{$hitB}{'length'} = $min_length;
					$best{$hitA}{$hitB}{'path'} = $path;
					#$best{$hitA}{'nodeB'} = $hitB;
				}
			}
			else {
				$best{$hitA}{$hitB}{'enrichment'} = $enrichment;
				$best{$hitA}{$hitB}{'score'} = $score;
				$best{$hitA}{$hitB}{'length'} = $min_length;
				$best{$hitA}{$hitB}{'path'} = $path;
				#$best{$hitA}{'nodeB'} = $hitB;
			}
			
		}
	}

	my %hit_counter = ();
	foreach my $a (keys %best) {
		foreach my $b (keys %{$best{$a}}) {
			if(exists $hits->{$a}) {
				$hit_counter{'hits'}{$a}++;
			}
			if(exists $hits->{$b}) {
				$hit_counter{'hits'}{$b}++;
			}
			foreach my $path_node (@{$best{$a}{$b}{'path'}}) {
				$hit_counter{'path'}{$path_node}++;
			}	
		}
	}
	print STDERR localtime() . "\tbest paths contain " . scalar(keys %{$hit_counter{'hits'}}) . 
		" hits and " . scalar(keys %{$hit_counter{'path'}}) . " path nodes\n";				
	
	return \%best;
}




sub evaluate_shortest_paths_original($$$) {
	my ($shortest, $hits, $enrichments) = @_;
	my %visits = ();
	my %best = ();
	# Compute how often a node is visited
	#print STDERR localtime() . "\t\tComputing node visit frequency using shortest paths\n";
	foreach my $hitA (sort keys %{$shortest->{'shortest'}}) {
		my @lengths = sort {$a <=> $b} keys %{$shortest->{'shortest'}{$hitA}};
		my $min_length = $lengths[0];
		my $start_neigh = 1;
		my $end_neigh = $min_length - 1;
		
		foreach my $hitB (sort keys %{$shortest->{'shortest'}{$hitA}{$min_length}}) {
			my $path = $shortest->{'shortest'}{$hitA}{$min_length}{$hitB} ;
			foreach my $node (@$path[$start_neigh..$end_neigh]) {
				$visits{$node}++;
			}
		}
	}
	# Now evaluate same length paths
	#print STDERR localtime() . "\t\tScoring paths\n";
	foreach my $hitA (sort keys %{$shortest->{'shortest'}}) {
		my @lengths = sort {$a <=> $b} keys %{$shortest->{'shortest'}{$hitA}};
		my $min_length = $lengths[0];
		my $start_neigh = 1;
		my $end_neigh = $min_length - 1;
		foreach my $hitB (sort keys %{$shortest->{'shortest'}{$hitA}{$min_length}}) {
			my $path = $shortest->{'shortest'}{$hitA}{$min_length}{$hitB};
			my $score = 0;
			my $enrichment = 1;
			foreach my $node (@$path[$start_neigh..$end_neigh]) {
				$score += 50 if(exists $hits->{$node});
				
				if(exists $visits{$node}) {
					$score += $visits{$node};
					
				}
				if(exists $enrichments->{$node}{'enrichment'}) {
					$enrichment *= $enrichments->{$node}{'enrichment'};
				}
				else {
					unless(exists $hits->{$node}) {
						$enrichment *= 0.001; 
					}
				}
			}
			if(exists $best{$hitA}{'enrichment'}) {
				my $prior_enrichment = $best{$hitA}{'enrichment'};
				if($enrichment > $prior_enrichment) {
					$best{$hitA}{'enrichment'} = $enrichment;
					$best{$hitA}{'score'} = $score;
					$best{$hitA}{'length'} = $min_length;
					$best{$hitA}{'path'} = $path;
					$best{$hitA}{'nodeB'} = $hitB;
				}
			}
			else {
				$best{$hitA}{'enrichment'} = $enrichment;
				$best{$hitA}{'score'} = $score;
				$best{$hitA}{'length'} = $min_length;
				$best{$hitA}{'path'} = $path;
				$best{$hitA}{'nodeB'} = $hitB;
			}			
		}	
	}
	my %hit_counter = ();
	foreach my $a (keys %best) {
		my $b = $best{$a}{'nodeB'};
		if(exists $hits->{$a}) {
			$hit_counter{$a}++;
		}
	}
	print STDERR localtime() . "\tbest paths contain " . scalar(keys %hit_counter) . " hits\n";				
	
	return \%best;
}

#
#	Remove non-hit nodes from the graph that are not articulation points
#	Articulation is recomputed in an iterative manner
#
sub remove_non_hit_non_articulation($$$) {
	my ($best, $graph, $hits) = @_;
	
	my %nodes = ();
	my %deleted = ();
	# Now add best path nodes:
	foreach my $hitA (sort keys %$best) {
		foreach my $hitB (sort keys %{$best->{$hitA}}) {
			$nodes{'any'}{$hitA}++;
			$nodes{'hit'}{$hitA}++;	
			$nodes{'any'}{$hitB}++;
			$nodes{'hit'}{$hitA}++;
		
			my $path = $best->{$hitA}{$hitB}{'path'};
			foreach my $path_node (@$path) {
				my $state = 'non-hit';
				if(exists $hits->{$path_node}) {
					$state = 'hit';
				}
				$nodes{$state}{$path_node}++;
				$nodes{'any'}{$path_node}++;
			}	
		}
	}
	my @vertices = sort keys %{$nodes{'any'}};
	
	my $hit_counter = 0;
	foreach my $a (@vertices) {
		if(exists $hits->{$a}) {
			$hit_counter++;
		}
	}
	print STDERR localtime() . "\tsubseting network with " . scalar(@vertices) . " nodes including $hit_counter hits\n";	
		
	# Now subset the graph
	my $subgraph = $graph->subgraph(\@vertices);
	$subgraph->biconnectivity_clear_cache;
	my $is_connected = $subgraph->is_connected;
	my @cc = $subgraph->connected_components();
	if($is_connected) {
		print STDERR localtime() . "\tsubgraph is connected\n";
	}
	else {
		print STDERR localtime() . "\tsubgraph components: " . scalar(@cc) . "\n";
	}
	
	
	my @articulation_points = $subgraph->articulation_points;
	foreach my $articulation_point (sort @articulation_points) {
		if(exists $nodes{'non-hit'}{$articulation_point}) {
			delete $nodes{'non-hit'}{$articulation_point};
			$nodes{'articulation'}{$articulation_point}++;
		}
	}
	my $node_counter = scalar(keys  %{$nodes{'non-hit'}});
	foreach my $node (sort keys %{$nodes{'non-hit'}}) {
		$node_counter--;
		$subgraph->biconnectivity_clear_cache;
		my @articulation_points = $subgraph->articulation_points;
		my $is_articulation;
		foreach my $ap (sort @articulation_points) {
			if($ap eq $node) {
				$is_articulation++;
			}
		}
		next if($is_articulation);
		my $copy_subgraph = $subgraph->copy_graph;
		$copy_subgraph->biconnectivity_clear_cache;
		print STDERR localtime() . "\tTesting deletion of non-articulation node $node_counter $node --> ";
		$copy_subgraph->delete_vertex($node);		
		my $is_connected = $copy_subgraph->is_connected;
		my @cc_before = $subgraph->connected_components();
		my @cc = $copy_subgraph->connected_components();
		if(scalar(@cc_before) == scalar(@cc)) {
			print STDERR " No change in graph components (" . scalar(@cc_before) .  ") --> ";
			$subgraph->delete_vertex($node);
			$deleted{$node}++;
			print STDERR " Deleting\n";
		}
		else {
			print STDERR " Increase in components: " . scalar(@cc) . " vs " .
				scalar(@cc_before) . " --> ";
			print STDERR " Cannot delete\n";
		}
		
		
		#$subgraph->delete_vertex($node);
		#$deleted{$node}++;
	}
	print STDERR localtime() . "\tDeleted " . scalar(keys %deleted) . " non-hit nodes: " .
								join(", ", (sort keys %deleted)) . "\n";
								
	my @nodes = $subgraph->vertices;
	$hit_counter = 0;
	foreach my $node (@nodes) {
		if(exists $hits->{$node}) {
			$hit_counter++;
		}
	}	
	my @edges = $subgraph->edges;
	print STDERR localtime() . "\tGenerated network with " . scalar(@nodes) . " nodes and " .
						scalar(@edges) . " edges containing $hit_counter hits\n";		

							
	return $subgraph;
	
}

sub remove_non_hit_non_articulation_original($$$) {
	my ($best, $graph, $hits) = @_;
	
	my %nodes = ();
	my %deleted = ();
	# Now add best path nodes:
	foreach my $hitA (sort keys %$best) {
		my $hitB = $best->{$hitA}{'nodeB'};
		$nodes{'any'}{$hitA}++;
		$nodes{'hit'}{$hitA}++;	
		$nodes{'any'}{$hitB}++;
		$nodes{'hit'}{$hitA}++;
		
		my $path = $best->{$hitA}{'path'};
		foreach my $path_node (@$path) {
			my $state = 'non-hit';
			if(exists $hits->{$path_node}) {
				$state = 'hit';
			}
			$nodes{$state}{$path_node}++;
			$nodes{'any'}{$path_node}++;
		}
	}
	my @vertices = sort keys %{$nodes{'any'}};
	
	my $hit_counter = 0;
	foreach my $a (@vertices) {
		if(exists $hits->{$a}) {
			$hit_counter++;
		}
	}
	print STDERR localtime() . "\tsubseting network with " . scalar(@vertices) . " nodes including $hit_counter hits\n";	
		
	# Now subset the graph
	my $subgraph = $graph->subgraph(\@vertices);
	$subgraph->biconnectivity_clear_cache;
	my $is_connected = $subgraph->is_connected;
	my @cc = $subgraph->connected_components();
	if($is_connected) {
		print STDERR localtime() . "\tsubgraph is connected\n";
	}
	else {
		print STDERR localtime() . "\tsubgraph components: " . scalar(@cc) . "\n";
	}
	
	
	my @articulation_points = $subgraph->articulation_points;
	foreach my $articulation_point (sort @articulation_points) {
		if(exists $nodes{'non-hit'}{$articulation_point}) {
			delete $nodes{'non-hit'}{$articulation_point};
			$nodes{'articulation'}{$articulation_point}++;
		}
	}
	foreach my $node (sort keys %{$nodes{'non-hit'}}) {
		$subgraph->biconnectivity_clear_cache;
		my @articulation_points = $subgraph->articulation_points;
		my $is_articulation;
		foreach my $ap (sort @articulation_points) {
			if($ap eq $node) {
				$is_articulation++;
			}
		}
		next if($is_articulation);
		my $copy_subgraph = $subgraph->copy_graph;
		$copy_subgraph->biconnectivity_clear_cache;
		print STDERR localtime() . "\tTesting deletion of non-articulation node $node\n";
		$copy_subgraph->delete_vertex($node);		
		my $is_connected = $copy_subgraph->is_connected;
		my @cc = $copy_subgraph->connected_components();
		if($is_connected) {
			print STDERR localtime() . "\tsubgraph copy is connected\n";
			$subgraph->delete_vertex($node);
			$deleted{$node}++;
			print STDERR localtime() . "\tDeleting non-articulation node $node\n";
		}
		else {
			print STDERR localtime() . "\tsubgraph copy components: " . scalar(@cc) . "\n";
			print STDERR localtime() . "\tCannot delete non-articulation node $node\n";
		}
		
		
		#$subgraph->delete_vertex($node);
		#$deleted{$node}++;
	}
	print STDERR localtime() . "\tDeleted " . scalar(keys %deleted) . " non-hit nodes: " .
								join(", ", (sort keys %deleted)) . "\n";
								
	my @nodes = $subgraph->vertices;
	$hit_counter = 0;
	foreach my $node (@nodes) {
		if(exists $hits->{$node}) {
			$hit_counter++;
		}
	}	
	my @edges = $subgraph->edges;
	print STDERR localtime() . "\tGenerated network with " . scalar(@nodes) . " nodes and " .
						scalar(@edges) . " edges containing $hit_counter hits\n";		

							
	return $subgraph;
	
}

#
#	Select approved, investigational, quality probe and 3D ligandable nodes that 
#	connect to at least two hits each
#
sub parse_druggability_file($$) {
	my ($best_druggability_file, $enrichments) = @_;
	
	# read in druggability file
	open(IN, $best_druggability_file) or die;
	my %druggability = ();
	while(defined(my $line = <IN>)) {
		chomp($line);
		my ($pac, $drug_class) = split(/\t/, $line);
		$druggability{$drug_class}{$pac}++;
	}
	close(IN);
	
	my %eligible_druggable = ();
	# First establish which are eligible (connect to two hits)
	foreach my $drug_node (keys %{$druggability{'Approved'}}) {
		next if(exists $hits->{$drug_node}); # no need to worry about these ones
		my $drug_hits = 0;
		if(exists $enrichments->{$drug_node}{'hits'}) {
			$drug_hits = scalar( keys %{$enrichments->{$drug_node}{'hits'}});
		} 
		next unless($drug_hits > 1); # need at least two neighbours
		foreach my $drug_hit (keys %{$enrichments->{$drug_node}{'hits'}}) {
			$eligible_druggable{'Approved'}{$drug_node}{$drug_hit}++;
		}	
	}
	print STDERR localtime() . "\t\tFound " . 
					scalar(keys %{$eligible_druggable{'Approved'}}) .
					" eligible Approved drug candidates\n";
		
	
	# Now check investigational (connect to two hits)
	foreach my $drug_node (keys %{$druggability{'Investigational'}}) {
		next if(exists $hits->{$drug_node}); # no need to worry about these ones
		my $drug_hits = 0;
		if(exists $enrichments->{$drug_node}{'hits'}) {
			$drug_hits = scalar( keys %{$enrichments->{$drug_node}{'hits'}});
		} 
		next unless($drug_hits > 1); # need at least two neighbours
		unless(exists $eligible_druggable{$drug_node}) {
			foreach my $drug_hit (keys %{$enrichments->{$drug_node}{'hits'}}) {
				$eligible_druggable{'Investigational'}{$drug_node}{$drug_hit}++;
			}	
		}	
	}
	print STDERR localtime() . "\t\tFound " . 
					scalar(keys %{$eligible_druggable{'Investigational'}}) .
					" eligible Investigational drug candidates\n";
	
	# Now check quality probes (connect to two hits)
	foreach my $drug_node (keys %{$druggability{'Quality_probe'}}) {
		next if(exists $hits->{$drug_node}); # no need to worry about these ones
		my $drug_hits = 0;
		if(exists $enrichments->{$drug_node}{'hits'}) {
			$drug_hits = scalar( keys %{$enrichments->{$drug_node}{'hits'}});
		} 
		next unless($drug_hits > 1); # need at least two neighbours
		unless(exists $eligible_druggable{$drug_node}) {
			foreach my $drug_hit (keys %{$enrichments->{$drug_node}{'hits'}}) {
				$eligible_druggable{'Quality_probe'}{$drug_node}{$drug_hit}++;
			}
		}	
	}
	print STDERR localtime() . "\t\tFound " . 
					scalar(keys %{$eligible_druggable{'Quality_probe'}}) .
					" eligible Quality probes candidates\n";

	# Now check 3D ligandable (connect to two hits)
	foreach my $drug_node (keys %{$druggability{'3D_ligandable'}}) {
		next if(exists $hits->{$drug_node}); # no need to worry about these ones
		my $drug_hits = 0;
		if(exists $enrichments->{$drug_node}{'hits'}) {
			$drug_hits = scalar( keys %{$enrichments->{$drug_node}{'hits'}});
		} 
		next unless($drug_hits > 1); # need at least two neighbours
		unless(exists $eligible_druggable{$drug_node}) {
			foreach my $drug_hit (keys %{$enrichments->{$drug_node}{'hits'}}) {
				$eligible_druggable{'3D_ligandable'}{$drug_node}{$drug_hit}++;
			}
		}	
	}
	print STDERR localtime() . "\t\tFound " . 
					scalar(keys %{$eligible_druggable{'3D_ligandable'}}) .
					" eligible 3D ligandable candidates\n";
					
	return \%eligible_druggable;
}

#
#	Add select druggable to the interactome as follows:
#	For approved add all
#	For investigational add all targeting hit pairs not already addressed 
#	For quality probe add all targeting hit pairs not already addressed 
#	For 3D ligandable add all targeting hit pairs not already addressed 
#	Generate and return subgraph from pre-existing hit network and new druggable nodes
#
sub add_druggablility_to_graph($$$) {
	my ($subgraph, $nr_stringent_interactome, $eligible_druggable) = @_;
	 
	my %nodes = ();
	
	my @subnet_nodes = $subgraph->vertices;	
	foreach my $subnet_node (sort @subnet_nodes) {
		$nodes{$subnet_node}++;
	}	
	my %seen_drug_hit_pairs = ();
	#Add eligible druggability nodes:
	foreach my $approved (sort keys %{$eligible_druggable->{'Approved'}}) {
		my @hits = sort keys  %{$eligible_druggable->{'Approved'}{$approved}};
		for(my $i = 0; $i < scalar(@hits); $i++) {
			for(my $j = $i; $j < scalar(@hits); $j++) {
				$seen_drug_hit_pairs{$hits[$i]}{$hits[$j]}{'Approved'}{$approved}++;
				$seen_drug_hit_pairs{$hits[$j]}{$hits[$i]}{'Approved'}{$approved}++;
				$nodes{$approved}++;
			}
		}		
	}
	# Add investigational unless better class interactions already exist
	foreach my $investigational (sort keys %{$eligible_druggable->{'Investigational'}}) {
		my @hits = sort keys  %{$eligible_druggable->{'Investigational'}{$investigational}};
		for(my $i = 0; $i < scalar(@hits); $i++) {
			for(my $j = $i; $j < scalar(@hits); $j++) {
				#next if(exists $seen_drug_hit_pairs{$hits[$i]}{$hits[$j]});
				#next if(exists $seen_drug_hit_pairs{$hits[$j]}{$hits[$i]});
				$seen_drug_hit_pairs{$hits[$i]}{$hits[$j]}{'Investigational'}{$investigational}++;
				$seen_drug_hit_pairs{$hits[$j]}{$hits[$i]}{'Investigational'}{$investigational}++;
				$nodes{$investigational}++;
			}
		}	
	}
	# Add Quality probes unless better class interactions already exist
	foreach my $qprobe (sort keys %{$eligible_druggable->{'Quality_probe'}}) {
		my @hits = sort keys  %{$eligible_druggable->{'Quality_probe'}{$qprobe}};
		for(my $i = 0; $i < scalar(@hits); $i++) {
			for(my $j = $i; $j < scalar(@hits); $j++) {
				next if(exists $seen_drug_hit_pairs{$hits[$i]}{$hits[$j]});
				next if(exists $seen_drug_hit_pairs{$hits[$j]}{$hits[$i]});
				$seen_drug_hit_pairs{$hits[$i]}{$hits[$j]}{'Quality_probe'}{$qprobe}++;
				$seen_drug_hit_pairs{$hits[$j]}{$hits[$i]}{'Quality_probe'}{$qprobe}++;
				$nodes{$qprobe}++;
			}
		}
	}
	# Add 3D-ligandable unless better class interactions already exist
	foreach my $ligandable (sort keys %{$eligible_druggable->{'3D_ligandable'}}) {
		my @hits = sort keys  %{$eligible_druggable->{'3D_ligandable'}{$ligandable}};
		for(my $i = 0; $i < scalar(@hits); $i++) {
			for(my $j = $i; $j < scalar(@hits); $j++) {
				next if(exists $seen_drug_hit_pairs{$hits[$i]}{$hits[$j]});
				next if(exists $seen_drug_hit_pairs{$hits[$j]}{$hits[$i]});
				$seen_drug_hit_pairs{$hits[$i]}{$hits[$j]}{'3D_ligandable'}{$ligandable}++;
				$seen_drug_hit_pairs{$hits[$j]}{$hits[$i]}{'3D_ligandable'}{$ligandable}++;
				$nodes{$ligandable}++;
			}
		}
	}
	# Now subset the graph
	my @vertices = sort keys %nodes;
	my $druggable_graph = $nr_stringent_interactome->subgraph(\@vertices);
	
	my @nodes = $druggable_graph->vertices;
	my $hit_counter = 0;
	foreach my $node (@nodes) {
		if(exists $hits->{$node}) {
			$hit_counter++;
		}
	}	
	my @edges = $druggable_graph->edges;
	print STDERR localtime() . "\tGenerated network with " . scalar(@nodes) . " nodes and " .
						scalar(@edges) . " edges containing $hit_counter hits\n";		

	return $druggable_graph;
}

#
#	Remove drug-drug edges that are not bridges
#	Then remove any other non-articulation, non-hit, non-drug nodes possible
#
sub compute_bridges_and_simplify($$$) {
	my ($graph, $eligible_druggable, $hits) = @_;
	$graph->biconnectivity_clear_cache;
	my @graph_edges = $graph->edges;
	my @graph_nodes = $graph->vertices;	
	my @bridges = $graph->bridges;
	# melt eligible into a drug node lookup
	my %drug_nodes = ();
	foreach my $class (keys %$eligible_druggable) {
		foreach my $drug_node (keys %{$eligible_druggable->{$class}}) {
			$drug_nodes{$drug_node}++;
		}	
	}
	my %deleted_drug_drug_edges = ();
	foreach my $graph_edge (@graph_edges) {
		my ($left, $right) = @$graph_edge;	
		my $status = 'not determined';
		if(exists $hits->{$left} || exists $hits->{$right}) {
			$status = 'hit edge';
		}
		else {
			foreach my $bridge (@bridges) {
				my ($b_left, $b_right) = @$bridge;
				if(($b_left eq $left) && ($b_right eq $right)) {
					$status = 'bridge edge';
				}
				elsif(($b_left eq $right) && ($b_right eq $left)) {
					$status = 'bridge edge';
				}
				#elsif($b_right eq $left) {
				#	$status = 'bridge edge';
				#}
				#elsif($b_right eq $right) {
				#	$status = 'bridge edge';
				#}
			}	
		} 
 		if(exists $drug_nodes{$left} && exists $drug_nodes{$right}) {
 			unless($status eq 'bridge edge' || $status eq 'hit edge') {
 				$status = 'drug - drug removable';
 				$graph->delete_edge($left, $right);
 				$graph->biconnectivity_clear_cache;
 				@bridges = $graph->bridges;
 				$deleted_drug_drug_edges{$left . " -- " . $right}++;
 				#print STDERR localtime() . "\t" . $left . "\t" . $right . "\t" . $status . " edge deleted\n";
 			}
 		}
	}
	print STDERR localtime() . "\tDeleted " . scalar(keys %deleted_drug_drug_edges) . 
					" drug - drug edges: " . join(", ", (sort keys %deleted_drug_drug_edges)) .
					"\n";
	# Find enriched path nodes
	my %enriched_path_nodes = ();
	foreach my $graph_node (@graph_nodes) {
		next if(exists $hits->{$graph_node});
		next if(exists $drug_nodes{$graph_node});
		$enriched_path_nodes{$graph_node}++;
	}
	$graph->biconnectivity_clear_cache;	
	my @articulation_points = $graph->articulation_points;
	foreach my $articulation_point (@articulation_points) {
		if(exists $enriched_path_nodes{$articulation_point}) {
			delete $enriched_path_nodes{$articulation_point};	
		}
	}
	my %deleted_non_articulation_node = ();
	foreach my $non_articulation_added (keys %enriched_path_nodes) {
		$graph->biconnectivity_clear_cache;
		my @articulation_points = $graph->articulation_points;
		my $is_articulation;
		foreach my $ap (@articulation_points) {
			if($ap eq $non_articulation_added) {
				$is_articulation++;
			}
		}
		next if($is_articulation);		
		$graph->delete_vertex($non_articulation_added);
		$deleted_non_articulation_node{$non_articulation_added}++;
	}
	print STDERR localtime() . "\tDeleted " . scalar(keys %deleted_non_articulation_node) .
		" non-articulation nodes: " . join(", ", (sort keys %deleted_non_articulation_node)) . "\n";
	my @nodes = $graph->vertices;
	my $hit_counter = 0;
	foreach my $node (@nodes) {
		if(exists $hits->{$node}) {
			$hit_counter++;
		}
	}	
	my @edges = $graph->edges;
	print STDERR localtime() . "\tGenerated network with " . scalar(@nodes) . " nodes and " .
						scalar(@edges) . " edges containing $hit_counter hits\n";		
	return $graph; 
}

#
#	Remove quality probe and 3D ligandable nodes that are not needed for connectivity
#
sub remove_qp_3d_not_needed($$$) {
	my ($graph, $eligible_druggable, $hits) = @_;

	my @graph_nodes = $graph->vertices;	
	$graph->biconnectivity_clear_cache;
	my @articulation_points = $graph->articulation_points;
	my %to_delete = ();
	my %deleted = ();
	foreach my $graph_node (@graph_nodes) {
		next if(exists $hits->{$graph_node});
		my $is_articulation;
		foreach my $ap (@articulation_points) {
			if($ap eq $graph_node) {
				$is_articulation++;
			}
		}
		next if($is_articulation);
		next if(exists $eligible_druggable->{'Approved'}{$graph_node});
		next if(exists $eligible_druggable->{'Investigational'}{$graph_node});	
		$to_delete{$graph_node}++;
		
	}
	# First remove challenging
	foreach my $to_delete_node (sort keys %to_delete) {
		next if(exists $eligible_druggable->{'Quality_probe'}{$to_delete_node});
		next if(exists $eligible_druggable->{'3D_ligandable'}{$to_delete_node});
		$graph->biconnectivity_clear_cache;
		my @articulation_points = $graph->articulation_points;
		my $is_articulation;
		foreach my $ap (@articulation_points) {
			if($ap eq $to_delete_node) {
				$is_articulation++;
			}
		}
		next if($is_articulation);
		$graph->delete_vertex($to_delete_node);
		$deleted{$to_delete_node}++;
		delete $to_delete{$to_delete_node};
	}
	# First remove 3D
	foreach my $to_delete_node (sort keys %to_delete) {
		next if(exists $eligible_druggable->{'Quality_probe'}{$to_delete_node});
		$graph->biconnectivity_clear_cache;
		my @articulation_points = $graph->articulation_points;
		my $is_articulation;
		foreach my $ap (@articulation_points) {
			if($ap eq $to_delete_node) {
				$is_articulation++;
			}
		}
		next if($is_articulation);
		$graph->delete_vertex($to_delete_node);
		$deleted{$to_delete_node}++;
		delete $to_delete{$to_delete_node};
	}
	# Then remove QP
	foreach my $to_delete_node (sort keys %to_delete) {
		$graph->biconnectivity_clear_cache;
		my @articulation_points = $graph->articulation_points;
		my $is_articulation;
		foreach my $ap (@articulation_points) {
			if($ap eq $to_delete_node) {
				$is_articulation++;
			}
		}
		next if($is_articulation);
		$graph->delete_vertex($to_delete_node);
		$deleted{$to_delete_node}++;
	}
	print STDERR localtime() . "\tDeleted " . scalar(keys %deleted) . 
		" non-articulation 3D/QP nodes: " . join(", ", (sort keys %deleted)) . "\n";
		
	my @nodes = $graph->vertices;
	my $hit_counter = 0;
	foreach my $node (@nodes) {
		if(exists $hits->{$node}) {
			$hit_counter++;
		}
	}	
	my @edges = $graph->edges;
	print STDERR localtime() . "\tGenerated network with " . scalar(@nodes) . " nodes and " .
						scalar(@edges) . " edges containing $hit_counter hits\n";	
	return $graph;
}

#
#	Generate SIF file from final network
#
sub generate_sif($$)  {
	my ($graph, $nr_interactions) = @_;
	my @graph_edges = $graph->edges;
	#foreach my $graph_edge (@graph_edges) {
	#	print STDERR join("\t",  @$graph_edge) . "\n";
	#}
	
	my @graph_nodes = $graph->vertices;
	#foreach my $graph_node (@graph_nodes) {
	#	print STDERR $graph_node . "\n";
	#}
	
	print STDERR localtime() ."\tGenerating sif with " . scalar(@graph_edges) . 
		" edges and " . scalar(@graph_nodes) . " nodes\n";
	my $sif_output = "";
	foreach my $graph_edge (@graph_edges) {
		my ($left, $right) = @$graph_edge;
		#print STDERR $left . "\t" . $right . "\t";
		if(exists $nr_interactions->{'OncoKB->'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'OncoKB->' . "\t" . $right . "\n";
			#print STDERR 'OncoKB->' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'OncoKB->'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'OncoKB->' . "\t" . $left . "\n";
			#print STDERR 'OncoKB<-' . "\n";
			next;
		}
	
		if(exists $nr_interactions->{'Signor->'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'Signor->' . "\t" . $right . "\n";
			#print STDERR 'Signor->' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'Signor->'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'Signor->' . "\t" . $left . "\n";
			#print STDERR 'Signor<-' . "\n";
			next;
		}
		
		if(exists $nr_interactions->{'Signor-|'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'Signor-|' . "\t" . $right . "\n";
			#print STDERR 'Signor-|' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'Signor-|'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'Signor-|' . "\t" . $left . "\n";
			#print STDERR 'Signor|-' . "\n";
			next;
		}

		if(exists $nr_interactions->{'reaction->'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'reaction->' . "\t" . $right . "\n";
			#print STDERR 'reaction->' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'reaction->'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'reaction->' . "\t" . $left . "\n";
			#print STDERR 'reaction<-' . "\n";
			next;
		}

		if(exists $nr_interactions->{'direct (LigInt)'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'direct (LigInt)' . "\t" . $right . "\n";
			#print STDERR 'direct->LIGINT' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'direct (LigInt)'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'direct (LigInt)' . "\t" . $left . "\n";
			#print STDERR 'direct<-LIGINT' . "\n";
			next;
		}
		
		if(exists $nr_interactions->{'direct (X-ray)'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'direct (X-ray)' . "\t" . $right . "\n";
			#print STDERR 'direct->XRAY' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'direct (X-ray)'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'direct (X-ray)' . "\t" . $left . "\n";
			#print STDERR 'direct<-XRAY' . "\n";
			next;
		}

		if(exists $nr_interactions->{'direct (nonPDB)'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'direct (nonPDB)' . "\t" . $right . "\n";
			#print STDERR 'direct->NONPDB' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'direct (nonPDB)'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'direct (nonPDB)' . "\t" . $left . "\n";
			#print STDERR 'direct<-NONPDB' . "\n";
			next;
		}
		
		if(exists $nr_interactions->{'transcriptional'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'transcriptional' . "\t" . $right . "\n";
			#print STDERR 'transcriptional->' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'transcriptional'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'transcriptional' . "\t" . $left . "\n";
			#print STDERR 'transcriptional<-' . "\n";
			next;
		}
		
		if(exists $nr_interactions->{'complex'}{$left}{$right}) {
			$sif_output .= $left . "\t" . 'complex' . "\t" . $right . "\n";
			#print STDERR 'complex->' . "\n";
			next;
		}
		elsif(exists $nr_interactions->{'complex'}{$right}{$left}) {
			$sif_output .= $right . "\t" . 'complex' . "\t" . $left . "\n";
			#print STDERR 'complex<-' . "\n";
			next;
		}
		print STDERR "\tNot found $left\t$right\n";	
	}	
	return $sif_output;
}
