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

# Get OS name
my $os_name = $^O;

# Turn on autoflush
$|++;

# Max number of forks to use
my $max_forks = get_free_cpus();

# Server port
my $port = 10003;

# Stores executing machine hostnames
my @machines;
my %machines;

# Path to text file containing computers to run on
my $machine_file_path;

# MrBayes block which will be used for each run
my $mb_block;

# Where this script is located 
my $script_path = abs_path($0);

# Directory script was called from
my $init_dir = abs_path(".");

# Where the script was called from
my $initial_directory = $ENV{PWD};

# General script settings
my $no_forks;

# Allow for reusing info from an old run
my $input_is_dir = 0;

# Allow user to specify mbsum output as input for the script
my $input_is_mbsum = 0;

# How the script was called
my $invocation = "perl bucky.pl @ARGV";

# Name of output directory
my $project_name = "bucky-".int(time());
#my $project_name = "bucky-dir";

# BUCKy settings
my $alpha = 1;
my $ngen = 1000000;

my @unlink;

# Read commandline settings
GetOptions(
	"no-forks"          => \$no_forks,
	"machine-file=s"    => \$machine_file_path,
	"alpha|a=s"         => \$alpha,
	"ngen|n=i"          => \$ngen,
	"port=i"            => \$port,
	"no-mbsum|s"        => \$input_is_mbsum,
	"n-threads|T=i"     => \$max_forks,
	"out-dir|o=s"       => \$project_name,
	"server-ip=s"       => \&client, # for internal usage only
	"help|h"            => sub { print &help; exit(0); },
	"usage"             => sub { print &usage; exit(0); },
);


# Get paths to required executables
my $bucky = check_path_for_exec("bucky");
my $mbsum = check_path_for_exec("mbsum") if (!$input_is_mbsum);

# Check that BUCKy version >= 1.4.4
check_bucky_version($bucky);

my $archive = shift(@ARGV);

# Some error checking
die "You must specify an archive file.\n\n", &usage if (!defined($archive));
die "Could not locate '$archive', perhaps you made a typo.\n" if (!-e $archive);
die "Could not locate '$machine_file_path'.\n" if (defined($machine_file_path) && !-e $machine_file_path);
die "Invalid alpha for BUCKy specified, input must be a float or 'infinity'.\n" if ($alpha !~ /(^inf(inity)?)|(^\d+(\.\d+)?$)/i);

# Input is a previous run directory, reuse information
$input_is_dir++ if (-d $archive);

# Determine which machines we will run the analyses on
if (defined($machine_file_path)) {

	# Get list of machines
	print "Fetching machine names listed in '$machine_file_path'...\n";
	open(my $machine_file, '<', $machine_file_path);
	chomp(@machines = <$machine_file>);
	close($machine_file);

	# Check that we can connect to specified machines
	foreach my $index (0 .. $#machines) {
		my $machine = $machines[$index];
		print "  Testing connection to: $machine...\n";

		# Attempt to ssh onto machine with a five second timeout
		my $ssh_test = `timeout 5 ssh -v $machine exit 2>&1`;

		# Look for machine's IP in test connection
		my $machine_ip;
		if ($ssh_test =~ /Connecting to \S+ \[(\S+)\] port \d+\./s) {
			$machine_ip = $1;
		}

		# Could connect but passwordless login not enabled
		if ($ssh_test =~ /Are you sure you want to continue connecting \(yes\/no\)/s) {
			print "    Connection to $machine failed, removing from list of useable machines (passwordless login not enabled).\n";
			splice(@machines, $index, 1);
		}
		# Successful connection
		elsif (defined($machine_ip)) {
			print "    Connection to $machine [$machine_ip] successful.\n";
			$machines{$machine} = $machine_ip;
		}
		# Unsuccessful connection
		else {
			print "    Connection to $machine failed, removing from list of useable machines.\n";
			splice(@machines, $index, 1);
		}
	}
}

print "\nScript was called as follows:\n$invocation\n\n";

my $archive_root;
my $archive_root_no_ext;
if (!$input_is_dir) {

	# Clean run with no prior output

	# Extract name information from input file

	if ($input_is_mbsum) {
		($archive_root = $archive) =~ s/.*\/(.*)/$1/;
		($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar(\.gz)?$)|(\.tgz$)/$2/;
	}
	else {
		($archive_root = $archive) =~ s/.*\/(.*)/$1/;
		($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.mb\.tar(\.gz)?$)|(\.mb\.tgz$)/$2/;
	}
	die "Could not determine archive root name, did you specify the proper input?\n" if (!defined($archive_root_no_ext));

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
#	($archive_root = $archive) =~ s/.*\/(.*)/$1/;
#	($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.mb\.tar(\.gz)?$)|(\.mb\.tgz$)/$2/;
	if ($input_is_mbsum) {
		($archive_root = $archive) =~ s/.*\/(.*)/$1/;
		($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar(\.gz)?$)|(\.tgz$)/$2/;
	}
	else {
		($archive_root = $archive) =~ s/.*\/(.*)/$1/;
		($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.mb\.tar(\.gz)?$)|(\.mb\.tgz$)/$2/;
	}
	die "Could not determine archive root name, is your input file properly named?\n" if (!defined($archive_root_no_ext));
}

# The name of the output archive
my $mbsum_archive = "$archive_root_no_ext.mbsum.tar.gz";
my $bucky_archive = "$archive_root_no_ext.BUCKy.tar";
my $quartet_output = "$archive_root_no_ext.CFs.csv";

$mbsum_archive = $archive if ($input_is_mbsum);

chdir($project_name);

# Change how Ctrl+C is interpreted to allow for clean up
$SIG{'INT'} = 'INT_handler';

# Define and initialize directories
my $mb_out_dir = "mb-out/";
my $mb_sum_dir = "mb-sum/";

mkdir($mb_out_dir) or die "Could not create '$mb_out_dir': $!.\n" if (!-e $mb_out_dir);
mkdir($mb_sum_dir) or die "Could not create '$mb_sum_dir': $!.\n" if (!-e $mb_sum_dir);

# Check if completed genes from a previous run exist
my %complete_quartets;
if (-e $bucky_archive && -e $quartet_output) {
	print "Archive containing completed quartets found for this dataset found in '$bucky_archive'.\n";
	print "Completed quartets within in this archive will be removed from the job queue.\n\n";

	# Because it takes longer to append to a tarball than append to a text file, the tarball most
	# likely has fewer quartet entries in it, we must therefore account for this ensuring that
	# the tarball and csv have the same quartet entries 

	# See which quartets in the tarball are complete
	chomp(my @complete_quartets_tarball = `tar tf '$init_dir/$project_name/$bucky_archive'`);

	# Add quartets to a hash for easier lookup
	my %complete_quartets_tarball;
	foreach my $complete_quartet (@complete_quartets_tarball) {
		$complete_quartets_tarball{$complete_quartet}++;
	}

	# Load csv into memory
	#open(my $quartet_output_file, "<", $quartet_output);
	open(my $quartet_output_file, "<", "$init_dir/$project_name/$quartet_output");
	(my @quartet_info = <$quartet_output_file>);
	close($quartet_output_file);

	# Remove header line
	my $header = shift(@quartet_info);

	# Rewrite csv to include only quartets also contained in tarball
	#open($quartet_output_file, ">", $quartet_output);
	open($quartet_output_file, ">", "$init_dir/$project_name/$quartet_output");
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

my @taxa;
my @genes;
if ($input_is_mbsum) {

	# Unarchive input genes 
	chomp(@genes = `tar xvf '$init_dir/$project_name/$archive' -C $mb_sum_dir 2>&1`);
	@genes = map { s/x //; $_ } @genes if ($os_name eq "darwin");

	die "No genes found in '$archive'.\n" if (!@genes);

	# Move into MrBayes output directory
	chdir($mb_sum_dir);

	# Remove subdirectories that may have been a part of the tarball
	my @dirs;
	my @gene_roots;
	foreach my $gene (@genes) {
		(my $gene_root = $gene) =~ s/.*\///;
		next if ($gene eq $gene_root);

		# Move mbsum files so they are no longer in subdirectories
		if (!-d $gene) {
			system("mv '$gene' '$gene_root'");
			push(@gene_roots, $gene_root);
		}
		else {
			push(@dirs, $gene);
			#splice(@genes, $index, 1);
		}
	}
	@genes = @gene_roots if (@gene_roots);

	# Remove any subdirectories that may have been in the input
	foreach my $dir (@dirs) {
		remove_tree($dir);
	}

	# Parse taxa present in each gene, determine which are shared across all genes
	my %taxa;
	foreach my $gene (@genes) {
		my @taxa = @{parse_mbsum_taxa($gene)};

		# Count taxa present
		foreach my $taxon (@taxa) {
			$taxa{$taxon}++;
		}
	}
	
	# Add taxa present in all genes to analysis
	foreach my $taxon (keys %taxa) {
		if ($taxa{$taxon} == scalar(@genes)) {
			push(@taxa, $taxon);
		}
	}

	# Archive genes
	system("tar", "czf", $mbsum_archive, @genes);
}
else {

	# Unarchive input genes 
	chomp(@genes = `tar xvf '$init_dir/$archive' -C '$mb_out_dir' 2>&1`);
	@genes = map { s/x //; $_ } @genes if ($os_name eq "darwin");

	# Move into MrBayes output directory
	chdir($mb_out_dir);

	# Check that each gene has a log file
	my %taxa;
	foreach my $gene (@genes) {

		# Unzip a single gene
		chomp(my @mb_files = `tar xvf '$gene' 2>&1`);
		@mb_files = map { s/x //; $_ } @mb_files if ($os_name eq "darwin");

		# Locate the log file output by MrBayes
		my $log_file_name;
		foreach my $file (@mb_files) {
			if ($file =~ /\.log$/) {
				$log_file_name = $file;
				last;
			}
		}
		die "Could not locate log file for '$gene'.\n" if (!defined($log_file_name));

		# Parse log file for run information
		my $mb_log = parse_mb_log($log_file_name);

		# Check for taxa present
		my @taxa = @{$mb_log->{TAXA}};
		foreach my $taxon (@taxa) {
			$taxa{$taxon}++;
		}

		# Clean up
		unlink(@mb_files);
	}

	# Add taxa present in all genes to analysis
	foreach my $taxon (keys %taxa) {
		if ($taxa{$taxon} == scalar(@genes)) {
			push(@taxa, $taxon);
		}
	}
}

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

print "Found ".scalar(@taxa)." taxa shared across all genes in this archive, ".scalar(@quartets).
	  " of $original_size possible quartets will be run using output from ".scalar(@genes)." total genes.\n";

# Go back to working directory
chdir("..");

# Determine whether or not we need to run mbsum on the specified input
my $should_summarize = 1;
if (-e $mbsum_archive && $input_is_dir && !$input_is_mbsum) {
	chomp(my @sums = `tar tf '$mbsum_archive'`) || die "Something appears to be wrong with '$mbsum_archive'.\n";

	# Check that each gene has actually been summarized, if not redo the summaries
	if (scalar(@sums) != scalar(@genes)) {
		unlink($mbsum_archive);
	}
	else {
		$should_summarize = 0;
	}
}

$should_summarize = 0 if ($input_is_mbsum);

# Summarize MrBayes output if needed
if ($should_summarize) {
	# Run mbsum on each gene
	print "Summarizing MrBayes output for ".scalar(@genes)." genes.\n";

	my @pids;
	foreach my $gene (@genes) {

		# Wait until a CPU is available
		until(okay_to_run(\@pids)) {};

		my $pid;
		until (defined($pid)) { $pid = fork(); usleep(30000); }

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
	#system("tar", "czf", $mbsum_archive, glob("$archive_root_no_ext*.sum"));
	system("tar", "czf", $mbsum_archive, glob("*.sum"));
	system("cp", $mbsum_archive, "..");
	chdir("..");
}

die "\nAll quartets have already been completed.\n\n" if (!@quartets);

# Returns the external IP address of this computer
chomp(my $server_ip = `dig +short myip.opendns.com \@resolver1.opendns.com 2>&1`);
if ($server_ip !~ /(?:[0-9]{1,3}\.){3}[0-9]{1,3}/) {
	print "Could not determine external IP address, only local clients will be created.\n";
	$server_ip = "127.0.0.1";
}

# Initialize a server
my $sock = IO::Socket::INET->new(
	LocalPort  => $port,
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
if (scalar(@machines) == 0) {
	push(@machines, $server_hostname);
	$machines{$server_hostname} = "127.0.0.1";
}
elsif (scalar(@machines) == 1) {
	# Check if the user input only the local machine in the config
	if ($machines{$machines[0]} eq $server_ip) {
		$machines{$machines[0]} = "127.0.0.1";
	}
}

my @pids;
foreach my $machine (@machines) {

	# Fork and create a client on the given machine
	my $pid = fork();	
	if ($pid == 0) {
		close(STDIN);
		close(STDOUT);
		close(STDERR);

		(my $script_name = $script_path) =~ s/.*\///;

		# Send over/copy files depending on where analyses will be run
		if ($machines{$machine} ne "127.0.0.1" && $machines{$machine} ne $server_ip) {

			# Send MrBayes summaries to remote machines
			if ($input_is_mbsum) {
				system("scp", "-q", "$mb_sum_dir/$mbsum_archive", $machine.":/tmp");
			}
			else {
				system("scp", "-q", $mbsum_archive, $machine.":/tmp");
			}

			# Send this script to the machine
			system("scp", "-q", $script_path, $machine.":/tmp");

			# Send BUCKy executable to the machine
			system("scp", "-q", $bucky, $machine.":/tmp");

			# Execute this perl script in client mode on the given machine
			# -tt forces pseudo-terminal allocation and lets us stop remote processes
			exec("ssh", "-tt", $machine, "perl", "/tmp/$script_name", $mbsum_archive, "--server-ip=$server_ip:$port");
		}
		else {
			# Send this script to the machine
			system("cp", $script_path, "/tmp");

			# Send BUCKy executable to the machine
			system("cp", $bucky, "/tmp");

			# Execute this perl script in client mode
			exec("perl", "/tmp/$script_name", "$init_dir/$project_name/$mb_sum_dir/$mbsum_archive", "--server-ip=127.0.0.1:$port");
		}

		exit(0);
	}
	else {
		push(@pids, $pid);
	}
}

# Move into mbsum directory
chdir($mb_sum_dir);

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
my %complete_queue;
my $complete_queue_max_size = 100;
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

		if (scalar(keys %complete_queue) > $complete_queue_max_size) {
			my $cwd = abs_path(".");
			chdir($mb_sum_dir);
			dump_quartets(\%complete_queue);
			chdir($cwd);
			undef(%complete_queue);
		}

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

				#push(@unlink, glob("$completed_quartet*"));
				push(@unlink, glob(abs_path(".")."/$completed_quartet*"));
				$complete_queue{$completed_quartet} = $quartet_statistics;
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
						print {$client} "CHDIR: ".abs_path(".")."\n";
					}

					# Invocation changes if we want to use a prior of infinity
					if ($alpha =~ /(^inf(inity)?)/i) {
						print {$client} "NEW: '$quartet' '--use-independence-prior -n $ngen'\n";
					}
					else {
						print {$client} "NEW: '$quartet' '-a $alpha -n $ngen'\n";
					}
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

# Dump remaining quartets
dump_quartets(\%complete_queue);

# Wait until all children have completed
foreach my $pid (@pids) {
	waitpid($pid, 0);
}

print "\n  All connections closed.\n";
print "Total execution time: ", sec2human(time() - $time), ".\n\n";

rmdir("$initial_directory/$project_name/$mb_sum_dir");

sub client {
	my ($opt_name, $address) = @_;	

	my ($server_ip, $port) = split(":", $address);

	chdir("/tmp");
	my $bucky = "/tmp/bucky";

	my $pgrp = $$;
	setpgrp();

	# Determine this host's IP
	chomp(my $ip = `dig +short myip.opendns.com \@resolver1.opendns.com`); 

	# Set IP to localhost if we don't have internet
	if ($ip !~ /(?:[0-9]{1,3}\.){3}[0-9]{1,3}/) {
		$ip = "127.0.0.1";
	}

	# Determine file name of mbsum archive the client should use
	my @ARGV = split(/\s+/, $invocation);
	shift(@ARGV); shift(@ARGV); # remove "perl" and "bucky.pl"
	my $mbsum_archive = shift(@ARGV);

	# Extract files from mbsum archive
	my @sums;
	if ($mbsum_archive =~ /\//) {
		chomp(@sums = `tar tf '$mbsum_archive'`);
	}
	else {
		chomp(@sums = `tar xvf '$mbsum_archive' 2>&1`);
		@sums = map { s/x //; $_ } @sums if ($os_name eq "darwin");
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
	$SIG{HUP} = sub { unlink($0, $bucky); unlink(@sums); unlink($mbsum_archive); kill -15, $$; };
	$SIG{TERM} = sub { unlink(glob("$quartet*")) if (defined($quartet)); exit(0); };

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

			# If client is local this needs to be defined now
			#chomp(@sums = `tar tf $mbsum_archive`) if (!@sums);

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
			system($bucky, split(" ", $bucky_settings), "-cf", 0, "-o", $quartet, "-p", $prune_file_path, @sums);
			unlink($prune_file_path);

			# Zip and tarball the results
			my @results = glob($quartet."*");
			my $quartet_archive_name = "$quartet.tar.gz";

			# Open concordance file and parse out the three possible resolutions
			my $num_genes = get_used_genes("$quartet.out");
			my $split_info = parse_concordance_output("$quartet.concordance", $num_genes);

			# Archive and compress results
			system("tar", "czf", $quartet_archive_name, @results);

			# Send the results back to the server if this is a remote client
			if ($server_ip ne "127.0.0.1" && $server_ip ne $ip) {
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
		#if ($server_ip ne $ip) {
		#	unlink(@sums, $mbsum_archive);
		#}
		#else {
		#	unlink(@sums);
		#}
	}

	exit(0);
}

sub dump_quartets {
	my $quartets = shift;

	my $pid;
	until (defined($pid)) { $pid = fork(); usleep(30000); }
	if ($pid == 0) {

		my @completed_quartets = keys %{$quartets};
		my @quartet_statistics = values %{$quartets};

		$SIG{TERM} = sub { close(STDIN); close(STDOUT); close(STDERR); unlink(@unlink); exit(0); };

		# Check if this is the first to complete, if so we must create CF output file
		if (!-e "../$quartet_output") {
			open(my $quartet_output_file, ">", "../$quartet_output");
			print {$quartet_output_file} "taxon1,taxon2,taxon3,taxon4,CF12.34,CF12.34_lo,CF12.34_hi,CF13.24,CF13.24_lo,CF13.24_hi,CF14.23,CF14.23_lo,CF14.23_hi,ngenes\n";
			foreach my $quartet (@quartet_statistics) {
				print {$quartet_output_file} $quartet,"\n";
			}
			close($quartet_output);
		}
		else {
			# Obtain a file lock on archive so another process doesn't simultaneously try to add to it
			open(my $quartet_output_file, ">>", "../$quartet_output");
			flock($quartet_output_file, LOCK_EX) || die "Could not lock '$quartet_output': $!.\n";
			seek($quartet_output_file, 0, SEEK_END) || die "Could not seek '$quartet_output': $!.\n";

			# Add completed gene
			#print {$quartet_output_file} @quartet_statistics,"\n";
			foreach my $quartet (@quartet_statistics) {
				print {$quartet_output_file} $quartet,"\n";
			}

			# Release lock
			flock($quartet_output_file, LOCK_UN) || die "Could not unlock '$quartet_output': $!.\n";
			close($quartet_output_file);
		}

		# Check if this is the first to complete, if so we must create BUCKy tarball
		if (!-e "../$bucky_archive") {
			system("tar", "cf", "../$bucky_archive", @completed_quartets);
			unlink(@completed_quartets);
		}
		else {
			my @unlink;
			foreach my $completed_quartet (@completed_quartets) {
				(my $quartet = $completed_quartet) =~ s/\.tar\.gz//;
				push(@unlink, glob("$quartet*"));
			}
			#$SIG{TERM} = sub { close(STDIN); close(STDOUT); close(STDERR); unlink(glob("$quartet*")); exit(0); };
			#$SIG{TERM} = sub { close(STDIN); close(STDOUT); close(STDERR); unlink(@unlink); exit(0); };

			# Obtain a file lock on archive so another process doesn't simultaneously try to add to it
			open(my $bucky_archive_file, "<", "../$bucky_archive");
			flock($bucky_archive_file, LOCK_EX) || die "Could not lock '$bucky_archive': $!.\n";

			# Add completed gene
			system("tar", "rf", "../$bucky_archive", @completed_quartets);
			unlink(@completed_quartets);

			# Release lock
			flock($bucky_archive_file, LOCK_UN) || die "Could not unlock '$bucky_archive': $!.\n";
			close($bucky_archive_file);
		}

		unlink(@unlink);
		undef(@unlink);

		exit(0);
	}
	else {
		push(@pids, $pid);
	}

#	unlink(@unlink);
#	undef(@unlink);

	return;
}

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
	$return .= ",$ngenes";

	return $return;
}

sub parse_mb_log {
	my $log_file_name = shift;

	# Open the specified mb log file and parse useful information from it

	my @taxa;
	my $ngen;
	my $nruns;
	my $burnin;
	my $burninfrac;
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
			$burninfrac = $1;
		}
		elsif ($line =~ /Setting chain burn-in to (\d+)/) {
			$burnin = $1;
		}
		elsif ($line =~ /Setting sample frequency to (\d+)/) {
			$samplefreq = $1;
		}
		elsif ($line =~ /Setting number of generations to (\d+)/) {
			$ngen = $1;
		}
	}
	close($log_file);

	if (defined($burnin)) {
		#return {'SAMPLEFREQ' => $samplefreq, 'NRUNS' => $nruns, 'BURNIN' => $burnin };
		return {'SAMPLEFREQ' => $samplefreq, 'NRUNS' => $nruns, 'BURNIN' => $burnin, 'TAXA' => \@taxa};
	}
	else {
		return {'NGEN' => $ngen, 'NRUNS' => $nruns, 'BURNINFRAC' => $burninfrac, 
				'SAMPLEFREQ' => $samplefreq, 'TAXA' => \@taxa};
	}
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

sub run_mbsum {
	my $tarball = shift;

	# Unzip specified tarball
	chomp(my @mb_files = `tar xvf '$mb_out_dir$tarball' -C '$mb_out_dir' 2>&1`);
	@mb_files = map { s/x //; $_ } @mb_files if ($os_name eq "darwin");

	# Determine name for this partition's log file
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

	# Determine number of trees mbsum should remove from each file

	my $trim;
	if ($mb->{BURNIN}) {
		$trim = $mb->{BURNIN} + 1;
	}
	else {
		$trim = ((($mb->{NGEN} / $mb->{SAMPLEFREQ}) * $mb->{NRUNS} * $mb->{BURNINFRAC}) / $mb->{NRUNS}) + 1;
	}

	# Summarize gene's tree files
	system("$mbsum '$mb_out_dir$gene_name.'*.t -n $trim -o '$mb_sum_dir$gene_name.sum' >/dev/null 2>&1");

	# Clean up extracted files
	chdir($mb_out_dir);
	unlink(@mb_files);
}

sub okay_to_run {
	my $pids = shift;

	# Free up a CPU by sleeping for 10 ms
	usleep(10000);

	my $current_forks = scalar(@{$pids});
	foreach my $index (reverse(0 .. $current_forks - 1)) {
		next if ($index < 0);

		my $pid = @{$pids}[$index];
		my $wait = waitpid($pid, WNOHANG);

		# Successfully reaped child
		if ($wait > 0) {
			$current_forks--;
			#splice(@pids, $index, 1);
			splice(@{$pids}, $index, 1);
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
	#dump_quartets(\%complete_queue);

	unlink(@unlink);
	
	# Kill ssh process(es) spawned by this script
	foreach my $pid (@pids) {
		#kill(-9, $pid);
		#kill(15, $pid);
		kill(1, $pid);
	}

	# Move into gene directory
	chdir("$initial_directory/$project_name");

	rmdir($mb_out_dir);

	# Try to delete directory once per second for five seconds, if it can't be deleted print an error message
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

sub sec2human {
	my $secs = shift;

	# Constants
	my $secs_in_min = 60;
	my $secs_in_hour = 60 * 60;
	my $secs_in_day = 24 * 60 * 60;

	$secs = int($secs);

	return "0 seconds" if (!$secs);

	# Calculate units of time
	my $days = int($secs / $secs_in_day);
	my $hours = ($secs / $secs_in_hour) % 24;
	my $mins = ($secs / $secs_in_min) % 60;
	$secs = $secs % 60;

	# Format return nicely
	my $time;
	if ($days) {
		$time .= ($days != 1) ? "$days days, " : "$days day, ";
	}
	if ($hours) {
		$time .= ($hours != 1) ? "$hours hours, " : "$hours hour, ";
	}
	if ($mins) {
		$time .= ($mins != 1) ? "$mins minutes, " : "$mins minute, ";
	}
	if ($secs) {
		$time .= ($secs != 1) ? "$secs seconds " : "$secs second ";
	}
	else {
		# Remove comma
		chop($time);
	}
	chop($time);

	return $time;
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
		chomp(@percent_free_cpu = `top -b -n2 -d0.05 | grep "Cpu(s)"`);
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
		$exec_path = abs_path($dir.$exec) if (-e $dir.$exec && -x $dir.$exec && !-d $dir.$exec);
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

sub check_bucky_version {
	my $bucky = shift;

	print "\nChecking for BUCKy version >= 1.4.4...\n";

	# Run BUCKy with --version and extract version info
	chomp(my @version_info = grep { /BUCKy version/ } `bucky --version`);
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

sub usage {
	return "Usage: bucky.pl [MRBAYES TARBALL]\n";
}

sub help {
print <<EOF; 
@{[usage()]}
Parallel execution of BUCKy on all possible quartets in a given alignment

  -a, --alpha            value of alpha to use when running BUCKy, use "infinity" for infinity (default: 1)      
  -n, --ngen             number of generations to run BUCKy MCMC chain (default: 1000000 generations)
  -o, --out-dir          name of the directory to store output files in (default: "bucky-" + Unix time of script invocation)
  -T, --n-threads        the number of forks ALL hosts running analyses can use concurrently (default: current number of free CPUs)
  -s, --no-mbsum         informs the script that the input is a tarball containing output already parsed by mbsum
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
