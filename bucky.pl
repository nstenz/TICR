#!/usr/bin/perl
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
use Time::HiRes qw(time usleep);

# Turn on autoflush
$|++;

# Max number of forks to use
my $max_forks = get_free_cpus();

# Server port
my $port = 10003;

# Stores executing machine hostnames
my @machines;

# Path to text file containing computers to run on
my $machine_file_path;

# MrBayes block which will be used for each run
my $mb_block;

# Where this script is located 
my $script_path = abs_path($0);

# Where the script was called from
my $initial_directory = $ENV{PWD};

# General script settings
my $no_forks;

# Allow for reusing info from an old run
my $input_is_dir = 0;

# How the script was called
my $invocation = "perl bucky.pl @ARGV";

# Name of output directory
my $project_name = "bucky-".time();
#my $project_name = "bucky-dir";

# BUCKy settings
my $alpha = 1;
my $ngen = 1000000;

# Read commandline settings
GetOptions(
	"no-forks"          => \$no_forks,
	"machine-file=s"    => \$machine_file_path,
	"alpha|a=f"         => \$alpha,
	"ngen|n=i"          => \$ngen,
	"port=i"            => \$port,
	"n-threads|T=i"     => \$max_forks,
	"out-dir|o=s"       => \$project_name,
	"server-ip=s"       => \&client, # for internal usage only
	"help|h"            => sub { print &help; exit(0); },
	"usage"             => sub { print &usage; exit(0); },
);


# Get paths to required executables
my $mbsum = check_path_for_exec("mbsum");
my $bucky = check_path_for_exec("bucky");

my $archive = shift(@ARGV);

# Some error checking
die "You must specify an archive file.\n\n", &usage if (!defined($archive));
die "Could not locate '$archive', perhaps you made a typo.\n" if (!-e $archive);
die "Could not locate '$machine_file_path'.\n" if (defined($machine_file_path) && !-e $machine_file_path);

# Input is a previous run directory, reuse information
$input_is_dir++ if (-d $archive);

# Determine which machines we will run the analyses on
if (defined($machine_file_path)) {
	print "Fetching machine names listed in '$machine_file_path'...\n";
	open(my $machine_file, '<', $machine_file_path);
	chomp(@machines = <$machine_file>);
	close($machine_file);
	print "  $_\n" foreach (@machines);
}

print "\nScript was called as follows:\n$invocation\n\n";

my $archive_root;
my $archive_root_no_ext;
if (!$input_is_dir) {

	# Clean run with no prior output

	# Extract name information from input file
	($archive_root = $archive) =~ s/.*\/(.*)/$1/;
	#($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar\.gz)|(\.tgz)/$2/;
	#($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar(\.gz)?$)|(\.tgz$)/$2/;
	#($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)((\.mb)?\.tar(\.gz)?$)|((\.mb)?\.tgz$)/$2/;
	#($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)((\.mb)?\.tar(\.gz)?$)|((\.mb)?\.tgz$)/$2/;
	($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.mb\.tar(\.gz)?$)|(\.mb\.tgz$)/$2/;

	# Initialize working directory
	# Remove conditional eventually
	mkdir($project_name) || die "Could not create '$project_name'$!.\n" if (!-e $project_name);

	my $archive_abs_path = abs_path($archive);
	# Remove conditional eventually
	run_cmd("ln -s $archive_abs_path $project_name/$archive_root") if (! -e "$project_name/$archive_root");
}
else {

	# Prior output available, set relevant variables

	$project_name = $archive;
	my @contents = glob("$project_name/*");

	# Determine the archive name by looking for a symlink
	my $found_name = 0;
	foreach my $file (@contents) {
		if (-l $file) {
			$file =~ s/\Q$project_name\E\///;
			$archive = $file;
			$found_name = 1;
		}
	}
	die "Could not locate archive in '$project_name'.\n" if (!$found_name);

	# Extract name information from input file
	($archive_root = $archive) =~ s/.*\/(.*)/$1/;
	#($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar\.gz)|(\.tgz)/$2/;
	#($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)((\.mb)?\.tar(\.gz)?$)|((\.mb)?\.tgz$)/$2/;
	($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.mb\.tar(\.gz)?$)|(\.mb\.tgz$)/$2/;
}

# The name of the output archive
my $mbsum_archive = "$archive_root_no_ext.mbsum.tar.gz";
my $bucky_archive = "$archive_root_no_ext.BUCKy.tar";
my $quartet_output = "$archive_root_no_ext.CFs.csv";

chdir($project_name);

# Change how Ctrl+C is interpreted to allow for clean up
$SIG{'INT'} = 'INT_handler';

# Define and initialize directories
my $mb_out_dir = "mb-out/";
my $mb_sum_dir = "mb-sum/";

mkdir($mb_out_dir) or die "Could not create '$mb_out_dir': $!.\n" if (!-e $mb_out_dir);
mkdir($mb_sum_dir) or die "Could not create '$mb_sum_dir': $!.\n" if (!-e $mb_sum_dir);

# Check if completed genes from a previous run exist
#my %complete_genes;
#if (-e $mb_archive) {
#	print "\nArchive containing completed MrBayes runs found for this dataset found in '$mb_archive'.\n";
#	print "Completed runs contained in this archive will be removed from the job queue.\n";
#
#	# Add gene names in tarball to list of completed genes
#	chomp(my @complete_genes = `tar tf $mb_archive`);
#	foreach my $gene (@complete_genes) {
#		$gene =~ s/\.tar\.gz//;
#		$complete_genes{$gene}++;
#	}
#}
my %complete_quartets;
if (-e $bucky_archive && -e $quartet_output) {
	print "\nArchive containing completed quartets found for this dataset found in '$bucky_archive'.\n";
	print "Completed quartets within in this archive will be removed from the job queue.\n";

	# Because it takes longer to append to a tarball than append to a text file, the tarball most
	# likely has fewer quartet entries in it, we must therefore account for this ensuring that
	# the tarball and csv have the same quartet entries 

	# See which quartets in the tarball are complete
	chomp(my @complete_quartets_tarball = `tar tf $bucky_archive`);

	# Add quartets to a hash for easier lookup
	my %complete_quartets_tarball;
	foreach my $complete_quartet (@complete_quartets_tarball) {
		$complete_quartets_tarball{$complete_quartet}++;
	}

	# Load csv into memory
	open(my $quartet_output_file, "<", $quartet_output);
	(my @quartet_info = <$quartet_output_file>);
	close($quartet_output_file);

	# Remove header line
	my $header = shift(@quartet_info);

	# Rewrite csv to include only quartets also contained in tarball
	open($quartet_output_file, ">", $quartet_output);
	print {$quartet_output_file} $header;
	foreach my $quartet (@quartet_info) {
		my ($taxon1, $taxon2, $taxon3, $taxon4) = split(",", $quartet);
		my $quartet_name = "$taxon1--$taxon2--$taxon3--$taxon4";
		my $quartet_name_tarball = $quartet_name.".tar.gz";

		if (exists($complete_quartets_tarball{$quartet_name_tarball})) {
			$complete_quartets{$quartet_name}++;
			print {$quartet_output_file} $quartet;
		}
	}
	close($quartet_output_file);
}

# Unarchive input genes 
chomp(my @genes = `tar xvf $archive -C $mb_out_dir`);

#chdir($gene_dir);
chdir($mb_out_dir);

#@genes = glob($archive_root_no_ext."*.nex");
#@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
#				(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
#				$a <=> $b } @genes;

# Unzip a single gene
chomp(my @mb_files = `tar xvf $genes[0]`);

# Locate the log file output by MrBayes
my $log_file_name;
foreach my $file (@mb_files) {
	if ($file =~ /\.log$/) {
		$log_file_name = $file;
		last;
	}
}
die "Could not locate log file for '$genes[0]'.\n" if (!defined($log_file_name));

# Parse log file for run information
my $mb_log = parse_mb_log($log_file_name);
my @taxa = @{$mb_log->{TAXA}};

# Create list of possible quartets
my @quartets = combine(\@taxa, 4);

my $original_size = scalar(@quartets);

# Remove completed quartets
if (%complete_quartets) {
	foreach my $index (reverse(0 .. $#quartets)) {
		my $quartet = $quartets[$index];
		$quartet = join("--", @{$quartet});
		if (exists($complete_quartets{$quartet})) {
			splice(@quartets, $index, 1);
		}
	}
}

print "Found ".scalar(@taxa)." taxa in this archive, ".scalar(@quartets)." of $original_size ".
      "possible quartets will be run using output from ".scalar(@genes)." total genes.\n";

# Go back to working directory
chdir("..");

# Determine whether or not we need to run mbsum on the specified input
my $should_summarize = 1;
if (-e $mbsum_archive && $input_is_dir) {
	chomp(my @sums = `tar tf $mbsum_archive`) || die "Something appears to be wrong with '$mbsum_archive'.\n";

	# Check that each gene has actually been summarized, if not redo the summaries
	if (scalar(@sums) != scalar(@genes)) {
		unlink($mbsum_archive);
	}
	else {
		$should_summarize = 0;
	}
}

if ($should_summarize) {
	# Run mbsum on each gene
	print "Summarizing MrBayes output for ".scalar(@genes)." genes.\n";

	my @pids;
	foreach my $gene (@genes) {

		# Wait until a CPU is available
		until(okay_to_run()) {};
		my $pid = fork();

		# The child fork
		if ($pid == 0) {
			#setpgrp();
			run_mbsum($gene);
			exit(0);
		}
		else {
			push(@pids, $pid);
		}
	}

	# Wait for all summaries to finish
	foreach my $pid (@pids) {
		waitpid($pid, 0);
	}
	undef(@pids);

	# Remove directory storing mb output
	remove_tree($mb_out_dir);

	# Archive and zip mb summaries
	chdir($mb_sum_dir);
	system("tar", "czf", $mbsum_archive, glob("$archive_root_no_ext*.sum"));
	system("cp", $mbsum_archive, "..");
	chdir("..");
}

die "\nAll quartets have already been completed.\n\n" if (!@quartets);

# Returns the external IP address of this computer
chomp(my $server_ip = `dig +short myip.opendns.com \@resolver1.opendns.com`);

# Initialize a server
my $sock = IO::Socket::INET->new(
	LocalAddr  => $server_ip.":".$port,
	Blocking   => 0,
	Reuse      => 1,
	Listen     => SOMAXCONN,
	Proto      => 'tcp') 
or die "Could not create server socket: $!.\n";
$sock->autoflush(1);

print "Job server successfully created.\n";

# Should probably do this earlier
# Determine server hostname and add to machines if none were specified by the user
chomp(my $server_hostname = `hostname`);
push(@machines, $server_hostname) if (scalar(@machines) == 0);

#@pids;
my @pids;
foreach my $machine (@machines) {

	# Fork and create a client on the given machine
	my $pid = fork();	
	if ($pid == 0) {
		close(STDIN);
		close(STDOUT);
		close(STDERR);

		(my $script_name = $script_path) =~ s/.*\///;

		# Send this script to the machine
		system("scp", "-q", $script_path, $machine.":/tmp");

		# Send BUCKy executable to the machine
		system("scp", "-q", $bucky, $machine.":/tmp");

		# Send MrBayes summaries to remote machines
		if ($machine ne $server_hostname) {
			system("scp", "-q", $mbsum_archive, $machine.":/tmp");
		}

		# Execute this perl script on the given machine
		# -tt forces pseudo-terminal allocation and lets us stop remote processes
		exec("ssh", "-tt", $machine, "perl", "/tmp/$script_name", $mbsum_archive, "--server-ip=$server_ip");
		exit(0);
	}
	else {
		push(@pids, $pid);
	}
}

# Move into mbsum directory
chdir($mb_sum_dir);

# Don't create zombies
#$SIG{CHLD} = 'IGNORE';

my $select = IO::Select->new($sock);

# Stores which job is next in queue 
my $job_number = 0;

# Number of open connections to a client
my $total_connections;

# Number of complete jobs (necessary?)
my $complete_count = 0;

# Number of connections server has closed
my $closed_connections = 0;

# Minimum number of connections server should expect
my $starting_connections = scalar(@machines);

my $time = time();
my $num_digits = get_num_digits({'NUMBER' => scalar(@quartets)});

# Begin the server's job distribution
while ((!defined($total_connections) || $closed_connections != $total_connections) || $total_connections < $starting_connections) {
	# Contains handles to clients which have sent information to the server
	my @clients = $select->can_read(0);

	# Free up CPU by sleeping for 10 ms
	usleep(10000);

	# Reap any children that we can
	foreach my $pid (@pids) {
		waitpid($pid, WNOHANG);
	}

	# Handle each ready client individually
	CLIENT: foreach my $client (@clients) {

		# Client requesting new connection
		if ($client == $sock) {
			$total_connections++;
			$select->add($sock->accept());
		}
		else {

			# Get client's message
			my $response = <$client>;
			next if (not defined($response)); # a response should never actually be undefined

			# Client wants to send us a file
			if ($response =~ /SEND_FILE: (.*)/) {
				my $file_name = $1;
				receive_file({'FILE_PATH' => $file_name, 'FILE_HANDLE' => $client});	
			}

			# Client has finished a job
			if ($response =~ /DONE '(.*)' '(.*)' \|\|/) {
				$complete_count++;				
				printf("  Analyses complete: %".$num_digits."d/%d.\r", $complete_count, scalar(@quartets));

				my $completed_quartet = $1;
				my $quartet_statistics = $2;

				# Check if this is the first to complete, if so we must create the archive
				if (!-e "../$bucky_archive") {
					system("tar", "cf", "../$bucky_archive", $completed_quartet, "--remove-files");
				}
				else {

					# Perform appending of new gene to tarball in a fork as this can take some time
					my $pid = fork();
					if ($pid == 0) {

						(my $quartet = $completed_quartet) =~ s/\.tar\.gz//;
						$SIG{TERM} = sub { close(STDIN); close(STDOUT); close(STDERR); unlink(glob("$quartet*")); exit(0); };

						# Obtain a file lock on archive so another process doesn't simultaneously try to add to it
						open(my $bucky_archive_file, "<", "../$bucky_archive");
						flock($bucky_archive_file, LOCK_EX) || die "Could not lock '$bucky_archive': $!.\n";

						# Add completed gene
						system("tar", "rf", "../$bucky_archive", $completed_quartet, "--remove-files");

						# Release lock
						flock($bucky_archive_file, LOCK_UN) || die "Could not unlock '$bucky_archive': $!.\n";
						close($bucky_archive_file);

						exit(0);
					}
					else {
						push(@pids, $pid);
					}
				}

				# Check if this is the first to complete, if so we must create CF output file
				if (!-e "../$quartet_output") {
					open(my $quartet_output_file, ">", "../$quartet_output");
					#print {$quartet_output_file} "taxon1\ttaxon2\ttaxon3\ttaxon4\tCF12|34\tCF13|24\tCF14|23\n";
					print {$quartet_output_file} "taxon1,taxon2,taxon3,taxon4,CF12.34,CF12.34_lo,CF12.34_hi,CF13.24,CF13.24_lo,CF13.24_hi,CF14.23,CF14.23_lo,CF14.23_hi\n";
					print {$quartet_output_file} $quartet_statistics,"\n";
					close($quartet_output);
				}
				else {

					# Perform appending of new quartet in a fork as this can take some time
					my $pid = fork();
					if ($pid == 0) {

						# Obtain a file lock on archive so another process doesn't simultaneously try to add to it
						open(my $quartet_output_file, ">>", "../$quartet_output");
						flock($quartet_output_file, LOCK_EX) || die "Could not lock '$quartet_output': $!.\n";
						seek($quartet_output_file, 0, SEEK_END) || die "Could not seek '$quartet_output': $!.\n";

						# Add completed gene
						print {$quartet_output_file} $quartet_statistics,"\n";

						# Release lock
						flock($quartet_output_file, LOCK_UN) || die "Could not unlock '$quartet_output': $!.\n";
						close($quartet_output_file);

						exit(0);
					}
					else {
						push(@pids, $pid);
					}
				}
			}

			# Client wants a new job
			if ($response =~ /NEW: (.*)/) {
				my $client_ip = $1;

				# Check if jobs remain in the queue
				if ($job_number < scalar(@quartets)) {
					printf("\n  Analyses complete: %".$num_digits."d/%d.\r", 0, scalar(@quartets)) if ($job_number == 0);

					my $quartet = join("--", @{$quartets[$job_number]});

					# Tell local clients to move into mbsum directory
					if ($client_ip eq $server_ip) {
						#print {$client} "CHDIR: ".abs_path($mb_sum_dir)."\n";
						print {$client} "CHDIR: ".abs_path(".")."\n";
					}
					print {$client} "NEW: '$quartet' '-a $alpha -n $ngen'\n";
					$job_number++;
				}
				else {
					# Client has asked for a job, but there are none remaining
					print {$client} "HANGUP\n";
					$select->remove($client);
					$client->close();
					$closed_connections++;
					next CLIENT;
				}
			}
		}
	}
}

# Wait until all children have completed
foreach my $pid (@pids) {
	waitpid($pid, 0);
}

print "\n  All connections closed.\n";
print "Total execution time: ", secs_to_readable({'TIME' => time() - $time}), "\n\n";

rmdir("$initial_directory/$project_name/$mb_sum_dir");

sub client {
	my ($opt_name, $server_ip) = @_;	

	chdir("/tmp");
	my $mb = "/tmp/mb";
	my $bucky = "/tmp/bucky";

	my $pgrp = $$;

	# Determine this host's IP
	chomp(my $ip = `dig +short myip.opendns.com \@resolver1.opendns.com`); 
	die "Could not establish an IP address for host.\n" if (not defined $ip);

	# Determine file name of mbsum archive the client should use
	my @ARGV = split(/\s+/, $invocation);
	shift(@ARGV); shift(@ARGV); # remove "perl" and "bucky.pl"
	my $mbsum_archive = shift(@ARGV);

	# Extract files from mbsum archive
	my @sums;
	if (-e $mbsum_archive) {
		chomp(@sums = `tar xvf $mbsum_archive`);
	}

	# Spawn more clients
	my @pids;
	my $total_forks = get_free_cpus(); 
	#my $total_forks = 1; 
	if ($total_forks > 1) {
		foreach my $fork (1 .. $total_forks - 1) {
			my $pid = fork();
			if ($pid == 0) {
				last;
			}
			else {
				push(@pids, $pid);
			}
		}
	}

	# The name of the quartet we are working on
	my $quartet;

	# Stores names of unneeded files
	my @unlink;

	# Change signal handling so killing the server kills these processes and cleans up
	#$SIG{CHLD} = 'IGNORE';
	#$SIG{HUP}  = sub { print "\n\n$$ received SIGHUP\n\n"; unlink($0, $bucky); unlink(@sums); unlink($mbsum_archive) if defined($mbsum_archive); kill -15, $$; };
	#$SIG{TERM} = sub { unlink(glob($quartet."*")) if defined($quartet); exit(0); };

	$SIG{HUP}  = sub {unlink($0, $bucky); unlink(@sums); unlink($mbsum_archive) if defined($mbsum_archive); kill -15, $$; };

	# Connect to the server
	my $sock = new IO::Socket::INET(
		PeerAddr  => $server_ip.":".$port,
		Proto => 'tcp') 
	or exit(0); 
	$sock->autoflush(1);

	print {$sock} "NEW: $ip\n";
	while (chomp(my $response = <$sock>)) {

		if ($response =~ /CHDIR: (.*)/) {
			chdir($1);
		}
		elsif ($response =~ /NEW: '(.*)' '(.*)'/) {
			$quartet = $1;
			my $bucky_settings = $2;

			#$SIG{TERM} = sub { close(STDIN); close(STDOUT); close(STDERR); unlink(@unlink); exit(0); };
			#$SIG{TERM} = sub { close(STDIN); close(STDOUT); close(STDERR); unlink(glob("$quartet*")); kill -9, $$; exit(0); };

			# If client is local this needs to be defined now
			chomp(@sums = `tar tf $mbsum_archive`) if (!@sums);
			#$SIG{TERM} = sub { unlink($0, $bucky); unlink(@sums); unlink($mbsum_archive); close(STDIN); close(STDOUT); close(STDERR); unlink(glob("$quartet*")); kill -9, $$; exit(0); };
			$SIG{TERM} = sub { close(STDIN); close(STDOUT); close(STDERR); unlink(glob("$quartet*")); kill -9, $$; exit(0); };

			# Create prune tree file contents required for BUCKy
			my $count = 0;
			my $prune_tree_output = "translate\n";
			foreach my $member (split("--", $quartet)) {
				$count++;
				$prune_tree_output .= " $count $member";
				if ($count == 4) {
					$prune_tree_output .= ";\n";
				}
				else {
					$prune_tree_output .= ",\n";
				}
			}

			# Write prune tree file
			my $prune_file_path = "$quartet-prune.txt";

			push(@unlink, $prune_file_path);

			open(my $prune_file, ">", $prune_file_path);
			print {$prune_file} $prune_tree_output;
			close($prune_file);

			push(@unlink, "$quartet.input", "$quartet.out", "$quartet.cluster", "$quartet.concordance", "$quartet.gene");

			# Run BUCKy on specified quartet
			system("$bucky $bucky_settings -cf 0 -o $quartet -p $prune_file_path @sums");
			unlink($prune_file_path);

			# Zip and tarball the results
			my @results = glob($quartet."*");
			my $quartet_archive_name = "$quartet.tar.gz";

			# Open concordance file and parse out the three possible resolutions
			my $split_info = parse_concordance_output("$quartet.concordance", scalar(@sums));

			# Archive and compress results
			#system("tar", "czf", $quartet_archive_name, @results, "--remove-files");
			system("tar", "czf", $quartet_archive_name, @results);

			# Send the results back to the server if this is a remote client
			if ($server_ip ne $ip) {
				send_file({'FILE_PATH' => $quartet_archive_name, 'FILE_HANDLE' => $sock});	
				unlink($quartet_archive_name);
			}

			unlink(@unlink);
			undef(@unlink);

			print {$sock} "DONE '$quartet_archive_name' '$split_info' || NEW: $ip\n";
		}
		elsif ($response eq "HANGUP") {
			last;
		}
	}

	# Have initial client wait for all others to finish and clean up
	if ($$ == $pgrp) {
		foreach my $pid (@pids) {
			waitpid($pid, 0);
		}
		unlink($0, $bucky);
		unlink(@sums, $mbsum_archive);
	}

	exit(0);
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

	return $return;
}

sub parse_mb_log {
	my $log_file_name = shift;

	# Open the specified mb log file and parse useful information from it

	my @taxa;
	my $ngen;
	my $nruns;
	my $burnin;
	my $samplefreq;
	open(my $log_file, "<", $log_file_name);
	while (my $line = <$log_file>) {
		if ($line =~ /Taxon\s+\d+ -> (\S+)/) {
			push(@taxa, $1);
		}
		elsif ($line =~ /Setting number of runs to (\d+)/) {
			$nruns = $1;
		}
		elsif ($line =~ /Setting burnin fraction to (\S+)/) {
			$burnin = $1;
		}
		elsif ($line =~ /Setting sample frequency to (\d+)/) {
			$samplefreq = $1;
		}
		elsif ($line =~ /Setting number of generations to (\d+)/) {
			$ngen = $1;
		}
		#last if ($line =~ /Exiting data block/);
	}
	close($log_file);

	return {'NGEN' => $ngen, 'NRUNS' => $nruns, 'BURNIN' => $burnin, 
	        'SAMPLEFREQ' => $samplefreq, 'TAXA' => \@taxa};
}

sub run_mbsum {
	my $tarball = shift;

	# Unzip specified tarball
	chomp(my @mb_files = `tar xvf $mb_out_dir$tarball -C $mb_out_dir`);

	my $log_file_name;
	foreach my $file (@mb_files) {
		if ($file =~ /\.log$/) {
			$log_file_name = $file;
			last;
		}
	}
	die "Could not locate log file for '$tarball'.\n" if (!defined($log_file_name));

	# Parse log file
	my $mb = parse_mb_log("$mb_out_dir$log_file_name");

	(my $gene_name = $tarball) =~ s/\.nex\.tar\.gz//;

	# Number of trees mbsum should remove from each file
	#my $trim = ((($ngen / $samplefreq) * $nruns * $burnin) / $nruns) + 1;
	my $trim = ((($mb->{NGEN} / $mb->{SAMPLEFREQ}) * $mb->{NRUNS} * $mb->{BURNIN}) / $mb->{NRUNS}) + 1;

	# Summarize gene's tree files
	system("$mbsum $mb_out_dir$gene_name.*.t -n $trim -o $mb_sum_dir$gene_name.sum >/dev/null 2>&1");
	#system("$mbsum $mb_out_dir$gene_name.*.t -n $trim -o $mb_sum_dir$gene_name.sum");

	# Clean up extracted files
	chdir($mb_out_dir);
	unlink(@mb_files);
}

sub okay_to_run {

	# Free up a CPU by sleeping for 10 ms
	usleep(10000);

	my $current_forks = scalar(@pids);
	foreach my $index (reverse(0 .. $#pids)) {
		my $pid = $pids[$index];
		my $wait = waitpid($pid, WNOHANG);

		# Successfully reaped child
		if ($wait > 0) {
			$current_forks--;
			splice(@pids, $index, 1);
		}
	}

	return ($current_forks < $max_forks);
}

sub hashsum {
	my $settings = shift;

	my $file_path = $settings->{'FILE_PATH'};

	open(my $file, "<", $file_path) or die "Couldn't open file '$file_path': $!.\n";
	my $md5 = Digest::MD5->new;
	my $md5sum = $md5->addfile(*$file)->hexdigest;
	close($file);

	return $md5sum;
}

sub send_file {
	my $settings = shift;

	my $file_path = $settings->{'FILE_PATH'};
	my $file_handle = $settings->{'FILE_HANDLE'};

	my $hash = hashsum({'FILE_PATH' => $file_path});
	print {$file_handle} "SEND_FILE: $file_path\n";

	open(my $file, "<", $file_path) or die "Couldn't open file '$file_path': $!.\n";
	while (<$file>) {
		print {$file_handle} $_;
	}
	close($file);

	print {$file_handle} " END_FILE: $hash\n";

	# Stall until we know status of file transfer
	while (defined(my $response = <$file_handle>)) {
		chomp($response);

		last if ($response eq "TRANSFER_SUCCESS");
		die "Unsuccessful file transfer, checksums did not match.\n" if ($response eq "TRANSFER_FAILURE");
	}
}

sub receive_file {
	my $settings = shift;

	my $file_path = $settings->{'FILE_PATH'};
	my $file_handle = $settings->{'FILE_HANDLE'};

	my $check_hash;
	open(my $file, ">", $file_path);
	while (<$file_handle>) {
		if ($_ =~ /(.*) END_FILE: (\S+)/) {
			print {$file} $1;
			$check_hash = $2;
			last;
		}
		else {
			print {$file} $_;
		}
	}
	close($file);

	# Use md5 hashsum to make sure transfer worked
	my $hash = hashsum({'FILE_PATH' => $file_path});
	if ($hash ne $check_hash) {
		die "Unsuccessful file transfer, checksums do not match.\n'$hash' - '$check_hash'\n"; # hopefully this never pops up
		print {$file_handle} "TRANSFER_FAILURE\n"
	}

	else {
		print {$file_handle} "TRANSFER_SUCCESS\n";
	}
}

sub INT_handler {

	# Kill ssh process(es) spawn by this script
	foreach my $pid (@pids) {
		#kill(-9, $pid);
		kill(15, $pid);
	}

	# Move into gene directory
	chdir("$initial_directory");

	# Try to delete directory five times, if it can't be deleted print an error message
	# I've found this method is necessary for analyses performed on AFS drives
	my $count = 0;
	until (!-e $mb_sum_dir || $count == 5) {
		$count++;

		remove_tree($mb_sum_dir, {error => \my $err});
		sleep(1);
	}
	#logger("Could not clean all files in './$gene_dir/'.") if ($count == 5);
	print "Could not clean all files in './$mb_sum_dir/'.\n" if ($count == 5);

	exit(0);
}

sub clean_up {
	my $settings = shift;

	my $remove_dirs = $settings->{'DIRS'};
	my $current_dir = getcwd();

#	chdir($alignment_root);
#	unlink(glob($gene_dir."$alignment_name*"));
#	#unlink($server_check_file) if (defined($server_check_file));
#
#	if ($remove_dirs) {
#		rmdir($gene_dir);
#	}
	chdir($current_dir);
}

sub get_num_digits {
	my $settings = shift;

	my $number = $settings->{'NUMBER'};

	my $digits = 1;
	while (floor($number / 10) != 0) {
		$number = floor($number / 10);
		$digits++;
	}

	return $digits;	
}

sub secs_to_readable {
	my $settings = shift;

	my %time;
	my $secs = $settings->{'TIME'};
	$time{'SEC'} = $secs;

	(my $decimal_secs = $secs) =~ s/.*(\.\d+)/$1/;

	my $mins = floor($secs / 60);
	if ($mins > 0) {
		$secs = $secs % 60;
		$secs += $decimal_secs if (defined $decimal_secs);

		$time{'SEC'} = $secs;
		$time{'MIN'} = $mins;

		my $hrs = floor($mins / 60);
		if ($hrs > 0) {
			$mins = $mins % 60;

			$time{'MIN'} = $mins;
			$time{'HOUR'} = $hrs;

			my $days = floor($hrs / 24);
			if ($days > 0) {
				$hrs  = $hrs % 24;	

				$time{'HOUR'} = $hrs;
				$time{'DAY'} = $days;
			}
		}
	}

	my $return;
	if (exists($time{'DAY'})) {
		$return = $time{'DAY'}." ".(($time{'DAY'} != 1) ? "days" : "day").
				  ", ".$time{'HOUR'}." ".(($time{'HOUR'} != 1) ? "hours" : "hour").
				  ", ".$time{'MIN'}." ".(($time{'MIN'} != 1) ? "minutes" : "minute").
				  ", ".$time{'SEC'}." ".(($time{'SEC'} != 1) ? "seconds" : "second").".";
	}
	elsif (exists($time{'HOUR'})) {
		$return = $time{'HOUR'}." ".(($time{'HOUR'} != 1) ? "hours" : "hour").
				  ", ".$time{'MIN'}." ".(($time{'MIN'} != 1) ? "minutes" : "minute").
				  ", ".$time{'SEC'}." ".(($time{'SEC'} != 1) ? "seconds" : "second").".";
	}
	elsif (exists($time{'MIN'})) {
		$return = $time{'MIN'}." ".(($time{'MIN'} != 1) ? "minutes" : "minute").
				  ", ".$time{'SEC'}." ".(($time{'SEC'} != 1) ? "seconds" : "second").".";
	}
	else {
		$return = $time{'SEC'}." ".(($time{'SEC'} != 1) ? "seconds" : "second").".";
	}
	return $return;
}

sub get_free_cpus {

	my $os_name = $^O;

	# Returns a two-member array containing CPU usage observed by top,
	# top is run twice as its first output is usually inaccurate
	my @percent_free_cpu;
	if ($os_name eq "darwin") {
		# Mac OS
		chomp(@percent_free_cpu = `top -i 1 -l 2 | grep "CPU usage"`);
	}
	else {
		# Linux
		chomp(@percent_free_cpu = `top -bn2d0.05 | grep "Cpu(s)"`);
	}

	my $percent_free_cpu = pop(@percent_free_cpu);

	if ($os_name eq "darwin") {
		# Mac OS
		$percent_free_cpu =~ s/.*?(\d+\.\d+)%\s+id.*/$1/;
	}
	else {
		# linux 
		$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
		$percent_free_cpu = eval($percent_free_cpu);
	}

	my $total_cpus;
	if ($os_name eq "darwin") {
		# Mac OS
		$total_cpus = `sysctl -n hw.ncpu`;
	}
	else {
		# linux
		$total_cpus = `grep --count 'cpu' /proc/stat` - 1;
	}

	my $free_cpus = ceil($total_cpus * $percent_free_cpu / 100);

	if ($free_cpus == 0 || $free_cpus !~ /^\d+$/) {
		$free_cpus = 1; # assume that at least one cpu can be used
	}
	
	return $free_cpus;
}

sub run_cmd {
	my $command = shift;

	my $return = system($command);

	if ($return) {
		logger("'$command' died with error: '$return'.\n");
		#kill(2, $parent_pid);
		exit(0);
	}
}

sub check_path_for_exec {
	my $exec = shift;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = abs_path($dir.$exec) if (-e $dir.$exec);
	}

	die "Could not find the following executable: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path));
	return $exec_path;
}

# I grabbed this from StackOverflow so that's why its style is different #DontFixWhatIsntBroken:
# https://stackoverflow.com/questions/10299961/in-perl-how-can-i-generate-all-possible-combinations-of-a-list
sub combine {
	my ($list, $n) = @_;
	die "Insufficient list members" if ($n > @$list);

	return map [$_], @$list if ($n <= 1);

	my @comb;

	for (my $i = 0; $i+$n <= @$list; ++$i) {
		my $val  = $list->[$i];
		my @rest = @$list[$i + 1 .. $#$list];
		push(@comb, [$val, @$_]) for combine(\@rest, $n - 1);
	}

	return @comb;
}

sub usage {
	return "Usage: bucky.pl [MRBAYES TARBALL]\n";
}

sub help {
print <<EOF; 
@{[usage()]}
Parallel execution of BUCKy on all possible quartets in a given alignment

  -a, --alpha            value of alpha to use when running BUCKy (default: 1)      
  -n, --ngen             number of generations to run BUCKy MCMC chain (default: 1000000 generations)
  -o, --out-dir          name of the directory to store output files in (default: "bucky-" + Unix time of script invocation)
  -T, --n-threads        the number of forks ALL hosts running analyses can use concurrently (default: current number of free CPUs)
  --machine-file         file name containing hosts to ssh onto and perform analyses on, passwordless login MUST be enabled
                         for each host specified in this file
  --port                 specifies the port to utilize on the server (Default: 10003)
  -h, --help             display this help and exit
  --usage                display proper script invocation format

Examples:
  perl bucky.pl align.mb.tar --machine-file hosts.txt     runs BUCKy using computers specified in hosts.txt using MrBayes output
                                                          stored in align.mb.tar

Mail bug reports and suggestions to <noah.stenz.github\@gmail.com>
EOF
exit(0);
}
