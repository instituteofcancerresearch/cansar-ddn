#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

sub read_hits_file($);
sub read_best_druggability_file($$);
sub read_sif_files($$$);
sub score_networks($);

my $dir;
my @sif_files;
my $output_file;
my $study;
my $best_druggability_file;
my $hits_file;
my $output_file;

GetOptions(
	"dir=s"		=>	\$dir,
	"sif=s"		=>	\@sif_files, 
	"hits=s"	=>	\$hits_file,
	"drug=s"	=>	\$best_druggability_file,
	"out=s"		=>	\$output_file,
);

			
my $data = read_hits_file($hits_file);
$data = read_best_druggability_file($data, $best_druggability_file);
$data = read_sif_files($data, $dir, \@sif_files);
$data = score_networks($data);
print "Top network: " . $data->{'network_best'} . "\n";

print "Network\tscore\tnodes\tedges\n";
foreach my $network (sort keys %{$data->{'report'}}) {
	print $network . "\t" . 
		$data->{'report'}{$network}{'node_score'} . "\t" .
		$data->{'report'}{$network}{'nodes'}  . "\t" .
		$data->{'report'}{$network}{'edges'} . "\n";
}		
exit;

sub score_networks($) {
	my $data = shift @_;
	foreach my $network (keys %{$data->{'network_input'}}) {
		my $node_score = 0;
		my $edges = 0;
		my %net_nodes = ();
		foreach my $hit (keys %{$data->{'network_input'}{$network}{'nodes'}{'hits'}}) {
			my $score = 1;
			my $status = $data->{'network_input'}{$network}{'nodes'}{'hits'}{$hit};
			if('Approved' eq $status) {
				$score = 10;
			}
			elsif('Investigational' eq $status) {
				$score = 8;
			}
			elsif('Quality_probe' eq $status) {
				$score = 4;
			}
			elsif('3D_ligandable' eq $status) {
				$score = 2;
			}
			$node_score += $score;
			$net_nodes{$hit}++;
		}
		foreach my $non_hit (keys %{$data->{'network_input'}{$network}{'nodes'}{'non_hits'}}) {
			my $score = 1;
			my $status = $data->{'network_input'}{$network}{'nodes'}{'non_hits'}{$non_hit};
			if('Approved' eq $status) {
				$score = 10;
			}
			elsif('Investigational' eq $status) {
				$score = 8;
			}
			elsif('Quality_probe' eq $status) {
				$score = 4;
			}
			elsif('3D_ligandable' eq $status) {
				$score = 2;
			}
			$node_score += $score;
			$net_nodes{$non_hit}++;
		}
		foreach my $node_a (keys %{$data->{'network_input'}{$network}{'edges'}}) {
			foreach my $node_b (keys %{$data->{'network_input'}{$network}{'edges'}{$node_a}}) {
				$edges++;
			}
		}
		$data->{'network_score'}{$node_score}{$network}{'nodes'} = scalar(keys %net_nodes);
		$data->{'network_score'}{$node_score}{$network}{'edges'} = $edges;
		$data->{'report'}{$network}{'node_score'} = $node_score;
		$data->{'report'}{$network}{'edges'} = $edges;
		$data->{'report'}{$network}{'nodes'} =  scalar(keys %net_nodes);
	}
	my @networks = sort {$b <=> $a} keys %{$data->{'network_score'}};
	my $top_score = $networks[0];
	my $best_net;
	my $min_node;
	my $min_edge;
	foreach my $net (keys %{$data->{'network_score'}{$top_score}}) {
		my $edges = $data->{'network_score'}{$top_score}{$net}{'edges'};
		my $nodes = $data->{'network_score'}{$top_score}{$net}{'nodes'};
		unless($best_net) {
			$best_net = $net;
			$min_node = $nodes;
			$min_edge = $edges;
			next;
		}
		if($nodes < $min_node) {
			$best_net = $net;
			$min_node = $nodes;
			$min_edge = $edges;
		}
		elsif($nodes == $min_node) {
			if($edges < $min_edge) {
				$best_net = $net;
				$min_node = $nodes;
				$min_edge = $edges;
			}
		}
	}
	$data->{'network_best'} = $best_net;
	#$data->{'network_best'}{$net}{'nodes'} = $min_node;
	#$data->{'network_best'}{$net}{'edges'} = $min_edge;
	return $data;
}

sub read_sif_files($$$) {
	my ($data, $dir, $files) = @_;
	foreach my $file (@$files) {
		open(IN, $dir . '/' . $file) or die;
		while(defined(my $line = <IN>)) {
			chomp($line);
			my ($a, $int, $b) = split(/\t/, $line);
			$data->{'network_input'}{$file}{'edges'}{$a}{$b} = $int;
			my $status_a = $data->{'best_druggability'}{$a};
			$status_a = 'Challenging' unless($status_a);
			if(exists $data->{'hits'}{$a}) {
				$data->{'network_input'}{$file}{'nodes'}{'hits'}{$a} = $status_a;
			}
			else {
				$data->{'network_input'}{$file}{'nodes'}{'non_hits'}{$a} = $status_a;
			}
			my $status_b = $data->{'best_druggability'}{$b};
			$status_b = 'Challenging' unless($status_b);
			if(exists $data->{'hits'}{$b}) {
				$data->{'network_input'}{$file}{'nodes'}{'hits'}{$b} = $status_b;
			}
			else {
				$data->{'network_input'}{$file}{'nodes'}{'non_hits'}{$b} = $status_b;
			}
		}
		close(IN);	
	}
	return $data;
}

sub read_best_druggability_file($$) {
	my ($data, $file) = @_;
	open(IN, $file) or die;
	while(defined(my $line = <IN>)) {
		chomp($line);
		my ($pac, $status) = split(/\t/, $line);
		$data->{'best_druggability'}{$pac} = $status;
		if(exists $data->{'hits'}{$pac}{'best_druggability'}) {
			$data->{'hits'}{$pac}{'best_druggability'} = $status;
		}
	}
	close(IN);
	return $data;
}

sub read_hits_file($) {
	my $file = shift @_;
	open(IN, $file) or die;
	my %data = ();
	while(defined(my $line = <IN>)) {
		chomp($line);
		$data{'hits'}{$line}{'best_druggability'} = 'Challenging';
	}
	close(IN);
	return \%data;
}


