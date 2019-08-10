#!/usr/bin/perl
## should be called:
## perl bucky-slurm.pl mbsum-folder -a alpha -n num-gen -q quartet-index -o output
## will use all files named mbsum-folder/*.in (one file = one gene)
##
## to test: perl bucky-slurm.pl mbsum -q 1
##
## Improvement over bucky.pl:
## - keeps all taxa (not only taxa shared in all genes)
## - only runs one quartet, so it can be used with any job scheduler
## things to do (fixit):
## - how many genes are acceptable per quartet?
## - what happens if a quartet has 0 genes? what does bucky do?
## - we do not want to do the translate table over and over, so this script
##   should check if it exists

use strict;
use warnings;
use POSIX;
use IO::Select;
use IO::Socket;
use Digest::MD5;
use Getopt::Long;
use Cwd qw(abs_path);
use Fcntl qw(:flock SEEK_END);
use POSIX qw(ceil :sys_wait_h);
use File::Path qw(remove_tree);
use File::Spec;
use Time::HiRes qw(time usleep);
##use integer;

# How the script was called
my $invocation = "perl bucky-slurm.pl @ARGV";
# Name of output directory
my $outdir = "bucky-".int(time());
# BUCKy settings
my $alpha = 1;
my $ngen = 1000000;
my $qind = 1;

# Read commandline settings
GetOptions(
	"alpha|a=s"         => \$alpha,
	"ngen|n=i"          => \$ngen,
	"quartet|q=i"          => \$qind,
	"out-dir|o=s"       => \$outdir,
	"help|h"            => sub { print &help; exit(0); },
	"usage"             => sub { print &usage; exit(0); },
);


## fixit: maybe we do not want to check this in every bucky, but how to do only once?
# Get paths to required executables
my $bucky = check_path_for_exec("bucky");
# Check that BUCKy version >= 1.4.4
check_bucky_version($bucky);

## this is the mbsum folder with files gene?.in
my $archive = shift(@ARGV);
# create output folder, if it doesn't exist already
mkdir($outdir) || die "Could not create '$outdir'$!.\n" if (!-e $outdir);

# Some error checking
die "You must specify an archive file.\n\n", &usage if (!defined($archive));
die "Could not locate '$archive', perhaps you made a typo.\n" if (!-e $archive);
die "Invalid alpha for BUCKy specified, input must be a float or 'infinity'.\n" if ($alpha !~ /(^inf(inity)?)|(^\d+(\.\d+)?$)/i);

print "\nScript was called as follows:\n$invocation\n\n";

opendir(BIN, $archive) or die "Can't open $archive: $!";
my @genes = grep { -T "$archive/$_" } readdir BIN;

# Parse taxa present in each gene
## fixit: create a file with all taxa: and then each script will check if this file exists, or not.
my %taxa;
foreach my $gene (@genes) {
    my @taxa = @{parse_mbsum_taxa("${archive}/$gene")};

    # Count taxa present
    foreach my $taxon (@taxa) {
	$taxa{$taxon}++;
    }
}



# Instead of keeping only taxa that are shared in all genes,
# we keep all taxa
my @taxa;
foreach my $taxon (keys %taxa) {
    #if ($taxa{$taxon} == scalar(@genes)) {
	push(@taxa, $taxon);
    #}
}

@taxa = sort @taxa;

my $ntax = scalar @taxa;
print "Read $ntax taxa\n";
my @quartet = whichQuartet($qind,$ntax);

# Check if concordance file exists
my $qscal_basename = join("-", @quartet);
my $qscal = File::Spec->catfile($outdir, $qscal_basename);
my $concordance_file = "$qscal.concordance";

## check also that the file is not empty
if(-e $concordance_file){
    print "BUCKy already run on this quartet: @quartet so nothing more to do\n";
    exit;
}

print "Calling BUCKy on the quartet: @quartet\n";

# Create prune tree file contents required for BUCKy
my $count = 0;
my $prune_tree_output = "translate\n";
foreach my $member (@quartet) {
    $count++;
    $prune_tree_output .= " $count $taxa[$member-1]";
    if ($count == 4) {
	$prune_tree_output .= ";\n";
    }
    else {
	$prune_tree_output .= ",\n";
    }
}

my $prune_file_path = "$qscal--prune.txt";
open(my $prune_file, ">", $prune_file_path);
print {$prune_file} $prune_tree_output;
close($prune_file);

#### Run BUCKy on specified quartet ####
## bucky is run in the *.in files of mbsum.
## CF cutoff for display          | -cf number                 | 0.05
## taxon set                      | -p prune-file              | common taxa
## output root file name          | -o name                    | run1
## skip genes with fewer taxa     | -sg                        | false

my $cmd = "$bucky -a $alpha -n $ngen -cf 0 -o $qscal -sg -p $prune_file_path $archive/*.in";
print "$cmd\n";
system($cmd);

#### Parsing output ####
# Open concordance file and parse out the three possible resolutions
my $num_genes = get_used_genes("$qscal.out");
my $keepQuartet = 1;
if($num_genes < 2){ ## fixit: how many genes is good enough?
    print "Quartet CF $qscal_basename were estimated with only $num_genes, will skip this quartet\n";
    $keepQuartet = 0;
}

if($keepQuartet){
    my $split_info = parse_concordance_output("$qscal.concordance", $num_genes);
    my $CFinfo_file = "$qscal.cf";
    open(my $CF_stream, ">", $CFinfo_file);
    print {$CF_stream} $split_info;
    close($CF_stream);
    print $split_info;
    print "\n";
}

# Delete extra output files
my @unlink = ("$qscal.cluster", "$qscal.gene", "$qscal.input", "$qscal.out", "$prune_file_path");
unlink(@unlink);
undef(@unlink);

sub get_used_genes {
	my $file_name = shift;

	# Number of genes actually used by BUCKy
	my $num_genes;

	# Open up the BUCKy output file
	open(my $bucky_out, "<", $file_name);
	while (my $line = <$bucky_out>) {
		if ($line =~ /Read (\d+) genes with a total of/) {
			$num_genes = $1;
			last;
		}
	}
	close($bucky_out);

	die "Error determining number of genes used by BUCKy ($file_name).\n" if (!defined($num_genes));

	return $num_genes;
}

sub parse_concordance_output {
	my ($file_name, $ngenes) = @_;

	my @taxa;
	my %splits;

	# Open up the specified output file
	open(my $concordance_file, "<", $file_name);

	my $split;
	my $in_translate;
	my $in_all_splits;
	while (my $line = <$concordance_file>) {

		# Parse the translate table
		if ($in_translate) {
			if ($line =~ /\d+ (.*)([,;])/) {
				my $taxon = $1;
				my $line_end = $2;
				push(@taxa, $taxon);

				$in_translate = 0 if ($line_end eq ';');
			}
		}

		# Parse the split information
		if ($in_all_splits) {

			# Set the split we are parsing information from
			if ($line =~ /^(\{\S+\})/) {
				my $current_split = $1;
				if ($current_split eq "{1,4|2,3}") {
					$split = "14|23";
				}
				elsif ($current_split eq "{1,3|2,4}") {
					$split = "13|24";
				}
				elsif ($current_split eq "{1,2|3,4}") {
					$split = "12|34";
				}
			}

			# Parse mean number of loci for split
			if ($line =~ /=\s+(\S+) \(number of loci\)/) {
				$splits{$split}->{"CF"} = $1 / $ngenes;
			}

			# Parse 95% confidence interval
			if ($line =~ /95% CI for CF = \((\d+),(\d+)\)/) {
				#$splits{$split}->{"95%_CI"} = "(".($1 / $ngenes).",".($2 / $ngenes).")";
				$splits{$split}->{"95%_CI_LO"} = ($1 / $ngenes);
				$splits{$split}->{"95%_CI_HI"} = ($2 / $ngenes);
			}
		}

		$in_translate++ if ($line =~ /^translate/);
		$in_all_splits++ if ($line =~ /^All Splits:/);
	}

	# Concat taxa names together
	my $return = join(",", @taxa);
	$return .= ",";

	# Concat split proportions with their 95% CI to return
	if (exists($splits{"12|34"})) {
		#$return .= $splits{"12|34"}->{"CF"}.$splits{"12|34"}->{"95%_CI"}."\t";
		$return .= $splits{"12|34"}->{"CF"}.",".$splits{"12|34"}->{"95%_CI_LO"}.",".$splits{"12|34"}->{"95%_CI_HI"}.",";
	}
	else {
		#$return .= "0(0,0)\t";
		$return .= "0,0,0,";
	}

	if (exists($splits{"13|24"})) {
		#$return .= $splits{"13|24"}->{"CF"}.$splits{"13|24"}->{"95%_CI"}."\t";
		$return .= $splits{"13|24"}->{"CF"}.",".$splits{"13|24"}->{"95%_CI_LO"}.",".$splits{"13|24"}->{"95%_CI_HI"}.",";
	}
	else {
		#$return .= "0(0,0)\t";
		$return .= "0,0,0,";
	}

	if (exists($splits{"14|23"})) {
		#$return .= $splits{"14|23"}->{"CF"}.$splits{"14|23"}->{"95%_CI"};
		$return .= $splits{"14|23"}->{"CF"}.",".$splits{"14|23"}->{"95%_CI_LO"}.",".$splits{"14|23"}->{"95%_CI_HI"};
	}
	else {
		#$return .= "0(0,0)";
		$return .= "0,0,0";
	}

	# Append number of genes used
	$return .= ",$ngenes\n";

	return $return;
}


sub parse_mbsum_taxa {
	my $mbsum_file_name = shift;

	# Open the specified mbsum file and parse its taxa list

	my @taxa;
	my $in_translate_block = 0;
	open(my $mbsum_file, "<", $mbsum_file_name);
	while (my $line = <$mbsum_file>) {
		$in_translate_block++ if ($line =~ /translate/);

		if ($in_translate_block == 1 && $line =~ /\d+\s+([^,;]+)/) {
		    push(@taxa, $1);
		}
	}
	close($mbsum_file);

	# Check if there were multiple translate blocks in the file which is indicative of an error with the creation of the file
	die "Something is amiss with '$mbsum_file_name', multiple translate blocks ($in_translate_block) were detected when there should only be one.\n" if ($in_translate_block > 1);

	# Check that we actually parsed something, otherwise input is improperly formatted/not mbsum output
	die "No taxa parsed for file '$mbsum_file_name', does '$archive' actually contain mbsum output?.\n" if (!@taxa);

	return \@taxa;
}

sub combination { ##http://stackoverflow.com/questions/10406411/efficient-inline-implementation-of-nck-in-perl
    my( $n, $r ) = @_;
    return unless defined $n && $n =~ /^\d+$/ && defined $r && $r =~ /^\d+$/;
    my $product = 1;
    while( $r > 0 ) {
        $product *= $n--;
        $product /= $r--;
    }
    my $rounded = sprintf "%.0f", $product;
    return $rounded;
}

sub whichQuartet {
    my $q = shift;
    my $n = shift;
    my $p = 4;
    die "  Error: number of taxa must be greater than 4.\n" if ($n<=4);
    my @quartet;
    while($n > 1)
    {
	#print "----\n";
	#print "n = $n \n";
	#print "q = $q \n";
	#print "p = $p \n";
        my $abs = combination($n-1,$p); #fixit: we don't want to compute this, we want to look for it in a table
	#print "abs = $abs \n";
	#my $subs = int($q-$abs);
	#if($subs > 0)
	if($q > $abs)
	{
	    #print "thinks q>abs \n";
	    push @quartet, $n;
	    $n = $n-1;
	    $p = $p-1;
	    $q = $q-$abs;
	}
	else
	{
	    $n = $n - 1;
	}
    }
    if(scalar @quartet == 3)
    {
	push @quartet, 1;
    }
    return sort @quartet;
}

sub usage {
	return "Usage: bucky.pl [MRBAYES TARBALL]\n";
}

sub help {
print <<EOF;
@{[usage()]}
Execution of BUCKy for one given quartet in a given alignment

  -a, --alpha            value of alpha to use when running BUCKy, use "infinity" for infinity (default: 1)
  -n, --ngen             number of generations to run BUCKy MCMC chain (default: 1000000 generations)
  -o, --out-dir          name of the directory to store output files in (default: "bucky-" + Unix time of script invocation)
  -q, --quartet          quartet index
  -h, --help             display this help and exit
  --usage                display proper script invocation format

Examples:
  perl bucky.pl mbsum-folder     runs BUCKy using mbsum output stored in the mbsum folder

Mail bug reports and suggestions to issues in github.com/TICR
EOF
exit(0);
}

sub check_bucky_version {
	my $bucky = shift;

	print "\nChecking for BUCKy version >= 1.4.4...\n";

	# Run BUCKy with --version and extract version info
	chomp(my @version_info = grep { /BUCKy version/ } `$bucky --version`);
	my $version_info = shift(@version_info);

	die "  Could not determine BUCKy version.\n" if (!defined($version_info));

	# Get the actual version number
	my $version;
	if ($version_info =~ /BUCKy\s+version\s+([^\s|,]+)/) {
		$version = $1;
	}
	die "  Could not determine BUCKy version.\n" if (!defined($version_info));

	# Version testing
	#my @versions = qw/2.1b 1.2.1000 1 0.9.8 2.3 1.4.5 1.4.3 1.4 1.500.2 1.4.4 1.4.4.1/;
	#foreach my $version (@versions) {

	print "  BUCKy version: $version.\n";

	# Split version number based on period delimiters
	my @version_parts = split(/\./, $version);

	# Future proofing if letters are ever used (we won't ever care about them)
	@version_parts = map { s/[a-zA-Z]+//g; $_ } @version_parts;
	die "  Error determining BUCKy version.\n" if (!@version_parts);

	# Check that version is >= 1.4.4
	if (defined($version_parts[0]) && $version_parts[0] > 1) {
		print "  BUCKy version check passed.\n";
		return;
	}
	elsif ((defined($version_parts[0]) && $version_parts[0] == 1) && (defined($version_parts[1]) && $version_parts[1] > 4)) {
		print "  BUCKy version check passed.\n";
		return;
	}
	elsif (((defined($version_parts[0]) && $version_parts[0] == 1) && (defined($version_parts[1]) && $version_parts[1] == 4))
		  && defined($version_parts[2]) && $version_parts[2] >= 4) {
		print "  BUCKy version check passed.\n";
		return;
	}
	else {
		die "  BUCKy version check failed, update to version >= 1.4.4.\n";
	}

	#print "\n";
	#}
}

sub check_path_for_exec {
	my $exec = shift;

	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = abs_path($dir.$exec) if (-e $dir.$exec && -x $dir.$exec && !-d $dir.$exec);
	}

	die "Could not find the following executable: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path));
	return $exec_path;
}

