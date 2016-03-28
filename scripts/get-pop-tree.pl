#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

# Use autoflush
$|++;

# Store script invocation
my $invocation = "perl get-pop-tree.pl @ARGV";

# Get QMC executable
my $qmc = get_QMC_exec();

# Load input
my $input = shift;

# Error checking
die "You must input a csv file containing quartet split information.\n" if (!defined($input));
die "Could not locate '$input'.\n" if (!-e $input);

# Print the current script settings
print "\nScript was called as follows:\n$invocation\n\n";

# Determine input root name
(my $input_root_no_ext = $input) =~ s/(.*\/)?(.*)\.CFs\.csv/$2/;

# File names for QMCN input and output
my $qmc_input = "$input_root_no_ext.QMC.txt";
my $qmc_output = "$input_root_no_ext.QMC.tre";

print "Parsing major resolution of each 4-taxon set... ";

# Parse input
my %taxa;
my @quartets;
open(my $input_file, "<", $input);
while (my $quartet = <$input_file>) {

	# Skip first line containing headers
	next if ($. == 1);

	# Split line
	my @line = split(",", $quartet);

	# Extract taxon names
	my ($taxon1, $taxon2, $taxon3, $taxon4) = ($line[0], $line[1], $line[2], $line[3]);
	$taxa{$taxon1}++; $taxa{$taxon2}++; $taxa{$taxon3}++; $taxa{$taxon4}++;

	# Extract split CFs
	my ($CF1234, $CF1324, $CF1423) = ($line[4], $line[7], $line[10]);

	# Put splits into an array to facilitate sorting
	my @splits = ($CF1234.":$taxon1,$taxon2|$taxon3,$taxon4",
	              $CF1324.":$taxon1,$taxon3|$taxon2,$taxon4",
	              $CF1423.":$taxon1,$taxon4|$taxon2,$taxon3");

	# Sort splits from highest to lowest CF
	my @sorted = sort { 
		(local $a = $a) =~ s/:.*//;	
		(local $b = $b) =~ s/:.*//;	
		$b <=> $a } @splits;

	# Extract sorted splits
	(my $CF_1 = $sorted[0]) =~ s/:.*//;
	(my $CF_2 = $sorted[1]) =~ s/:.*//;
	(my $CF_3 = $sorted[2]) =~ s/:.*//;
	(my $split_1 = $sorted[0]) =~ s/.*://;
	(my $split_2 = $sorted[1]) =~ s/.*://;
	(my $split_3 = $sorted[2]) =~ s/.*://;

	# All three quartets are equal, include all
	if ($CF_1 == $CF_2 && $CF_2 == $CF_3) {
		push(@quartets, $split_1);	
		push(@quartets, $split_2);	
		push(@quartets, $split_3);	
	}
	# Two quartets are equal, include both
	elsif ($CF_1 == $CF_2) {
		push(@quartets, $split_1);	
		push(@quartets, $split_2);	
	}
	# Only include the highest quartet
	else {
		push(@quartets, $split_1);	
	}
}
close($input_file);

print "done.\n";

# Create hash to allow for renaming of quartets based on sorted taxa names
my $id = 1;
my %taxon_to_id;
foreach my $taxon (sort {$a cmp $b} keys %taxa) {
	$taxon_to_id{$taxon} = $id;
	$id++;
}

# Rewrite quartets to include ids, not actual taxa names
foreach my $quartet (@quartets) {
	if ($quartet =~ /(\S+),(\S+)\|(\S+),(\S+)/) {
		$quartet = "$taxon_to_id{$1},$taxon_to_id{$2}|$taxon_to_id{$3},$taxon_to_id{$4}";
	}
}

# Output quartets to file
open(my $qmc_input_file, ">", $qmc_input);
print {$qmc_input_file} join(" ", @quartets),"\n";
close($qmc_input_file);

print "Running Quartet Max Cut...\n";

# Run Quartet Max Cut
system($qmc, "qrtt=$qmc_input", "otre=$qmc_output");
unlink($qmc_input);

# Open Quartet Max Cut output and replace taxa ids with actual names
open(my $qmc_output_file, "<", $qmc_output);
chomp(my $tree = <$qmc_output_file>);
close($qmc_output_file);

my %id_to_taxon = reverse(%taxon_to_id);
$tree =~ s/(\d+)/$id_to_taxon{$1}/g;

# Reopen Quartet Max Cut output and replace it with the modified tree
open($qmc_output_file, ">", $qmc_output);
print {$qmc_output_file} $tree,"\n";
close($qmc_output_file);

print "Quartet Max Cut complete, tree located in '$qmc_output'.\n\n";

sub get_QMC_exec {
	my $OS = $^O;
	my %execs = ('linux' => 'find-cut-Linux-64', 'darwin' => 'find-cut-Mac');

	my $exec = $execs{$OS};

	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = abs_path($dir.$exec) if (-e $dir.$exec);
	}

	die "Could not find the following executable to run Quartet MaxCut: '$exec'. This script requires this program in your path.\nQuartet MaxCut can be downloaded here: http://research.haifa.ac.il/~ssagi/software/QMCN.tar.gz" if (!defined($exec_path));
	return $exec_path;
}
