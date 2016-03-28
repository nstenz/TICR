#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use IO::Select;
use IO::Socket;
use Digest::MD5;
use Getopt::Long;
use Cwd qw(abs_path);
use Fcntl qw(:flock);
use File::Path qw(remove_tree);
use Time::HiRes qw(time usleep);

my $os_name = $^O;

# Turn on autoflush
$|++;

# Maximum number of threads to use
my $max_forks;

# Server port
my $port = 10002;

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

# Allow for reusing info from an old run
my $input_is_dir = 0;

# How the script was called
my $invocation = "perl mb.pl @ARGV";

# Name of output directory
my $project_name = "mb-".int(time());
#my $project_name = "mb-dir";

# Read commandline settings
GetOptions(
	#"no-forks"          => \$no_forks,
	"mb-block|m:s" => \$mb_block,
	"machine-file:s"    => \$machine_file_path,
	"check|c:f"         => \&check_nonconvergent,
	"remove|r:f"        => \&remove_nonconvergent,
    "out-dir|o=s"       => \$project_name,
	"n-threads|T"       => \$max_forks,
	"port=i"            => \$port,
	"server-ip:s"       => \&client, # for internal usage only
	"help|h"            => sub { print &help; exit(0); },
	"usage"             => sub { print &usage; exit(0); },
);


# Get paths to required executables
my $mb = check_path_for_exec("mb");

my $archive = shift(@ARGV);

# Some error checking
die "You must specify an archive file.\n\n", &usage if (!defined($archive));
die "Could not locate '$archive', perhaps you made a typo.\n" if (!-e $archive);
die "You specified a MrBayes run archive instead of an MDL gene archive.\n" if ($archive =~ /\.mb\.tar$/);
die "Could not locate '$machine_file_path'.\n" if (defined($machine_file_path) && !-e $machine_file_path);
die "You must specify a file containing a valid MrBayes block which will be appended to each gene.\n\n", &usage if (!defined($mb_block));
die "Could not locate '$mb_block', perhaps you made a typo.\n\n" if (!-e $mb_block);

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

print "\nScript was called as follows:\n$invocation\n";

# Load MrBayes block into memory
open(my $mb_block_file, "<", $mb_block) or die "Could not open '$mb_block': $!.\n";
my @mb_block = <$mb_block_file>;
close($mb_block_file);

my $archive_root;
my $archive_root_no_ext;
if (!$input_is_dir) {

	# Clean run with no prior output

	# Extract name information from input file
	($archive_root = $archive) =~ s/.*\/(.*)/$1/;
	($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar\.gz)|(\.tgz)/$2/;

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
			#$archive = $file;
			$archive = "$project_name/$file";
			$found_name = 1;
		}
	}
	die "Could not locate archive in '$project_name'.\n" if (!$found_name);

	# Extract name information from input file
	($archive_root = $archive) =~ s/.*\/(.*)/$1/;
	($archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar\.gz)|(\.tgz)/$2/;
}

# The name of the output archive
my $mb_archive = "$archive_root_no_ext.mb.tar";

chdir($project_name);

# Change how Ctrl+C is interpreted to allow for clean up
$SIG{'INT'} = 'INT_handler';

# Define and initialize directories
my $gene_dir = "genes/";
mkdir($gene_dir) or die "Could not create '$gene_dir': $!.\n" if (!-e $gene_dir);

# Check if completed genes from a previous run exist
my %complete_genes;
if (-e $mb_archive) {
	print "\nArchive containing completed MrBayes runs found for this dataset found in '$mb_archive'.\n";
	print "Completed runs contained in this archive will be removed from the job queue.\n";

	# Add gene names in tarball to list of completed genes
	chomp(my @complete_genes = `tar tf '$mb_archive'`);
	foreach my $gene (@complete_genes) {
		$gene =~ s/\.tar\.gz//;
		$complete_genes{$gene}++;
	}
}

# Unarchive input genes 
chomp(my @genes = `tar xvf '$init_dir/$archive' -C $gene_dir 2>&1`);
@genes = map { s/x //; $_ } @genes if ($os_name eq "darwin");

chdir($gene_dir);

# Remove completed genes
if (%complete_genes) {
	foreach my $index (reverse(0 .. $#genes)) {
		if (exists($complete_genes{$genes[$index]})) {
			unlink($genes[$index]);
			splice(@genes, $index, 1);
		}
	}
}

die "\nAll jobs have already completed.\n\n" if (!@genes);

# Append given MrBayes block to the end of each gene
print "\nAppending MrBayes block to each gene... ";
foreach my $gene (@genes) {
	open(my $gene_file, ">>", $gene) or die "Could not open '$gene': $!.\n";
	print {$gene_file} "\n", @mb_block;
	close($gene_file);
}
print "done.\n\n";

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

		# Move required datafiles to machines, initialize clients
		if ($machines{$machine} ne "127.0.0.1" && $machines{$machine} ne $server_ip) {
			# Send this script to the machine
			system("scp", "-q", $script_path, $machine.":/tmp");

			# Send MrBayes executable to the machine
			system("scp", "-q", $mb, $machine.":/tmp");

			# Execute this perl script on the given machine
			# -tt forces pseudo-terminal allocation and lets us stop remote processes
			exec("ssh", "-tt", "$machine", "perl", "/tmp/$script_name", "--server-ip=$server_ip:$port");
		}
		else {
			# Send this script to the machine
			system("cp", $script_path, "/tmp");

			# Send MrBayes executable to the machine
			system("cp", $mb, "/tmp");

			# Execute this perl script on the given machine
			exec("perl", "/tmp/$script_name", "--server-ip=127.0.0.1:$port");
		}

		exit(0);
	}
	else {
		push(@pids, $pid);
	}
}

#chdir($gene_dir);

my $select = IO::Select->new($sock);

# Don't create zombies
$SIG{CHLD} = 'IGNORE';

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
my $num_digits = get_num_digits({'NUMBER' => scalar(@genes)});

# Begin the server's job distribution
while ((!defined($total_connections) || $closed_connections != $total_connections) || $total_connections < $starting_connections) {
	# Contains handles to clients which have sent information to the server
	my @clients = $select->can_read(0);

	# Free up CPU by sleeping for 10 ms
	usleep(10000);

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
			if ($response =~ /DONE (.*) \|\|/) {
				$complete_count++;				
				printf("  Analyses complete: %".$num_digits."d/%d.\r", $complete_count, scalar(@genes));

				# Perform appending of new gene to tarball in a fork as this can take some time
				my $pid;
				until (defined($pid)) { $pid = fork(); usleep(30000); }

				if ($pid == 0) {

					# Check if this is the first to complete, if so we must create the directory
					my $completed_gene = $1;
					if (!-e "../$mb_archive") {
						system("touch", "$mb_archive");
						system("tar", "cf", "../$mb_archive", $completed_gene);
						unlink($completed_gene);
					}
					else {

						# Obtain a file lock on archive so another process doesn't simultaneously try to add to it
						open(my $mb_archive_file, "<", "../$mb_archive");
						flock($mb_archive_file, LOCK_EX) || die "Could not lock '$mb_archive_file': $!.\n";

						# Add completed gene
						system("tar", "rf", "../$mb_archive", $completed_gene);
						unlink($completed_gene);

						# Release lock
						flock($mb_archive_file, LOCK_UN) || die "Could not unlock '$mb_archive_file': $!.\n";
						close($mb_archive_file);

					}
					exit(0);
				}
				else {
					push(@pids, $pid);
				}
			}

			# Client wants a new job
			if ($response =~ /NEW: (.*)/) {
				my $client_ip = $1;

				# Check if jobs remain in the queue
				if ($job_number < scalar(@genes)) {
					printf("\n  Analyses complete: %".$num_digits."d/%d.\r", 0, scalar(@genes)) if ($job_number == 0);

					my $gene = $genes[$job_number];

					# Check whether the client is remote or local, send it needed files if remote
					if ($client_ip ne $server_ip) {

						# Fork to perform the file transfer and prevent stalling the server
						my $pid;
						until (defined($pid)) { $pid = fork(); usleep(30000); }

						#my $pid = fork(); 
						if ($pid == 0) {
							send_file({'FILE_PATH' => $gene, 'FILE_HANDLE' => $client});	
							unlink($gene);

							print {$client} "NEW: $gene\n";
							exit(0);
						}
						else {
							push(@pids, $pid);
						}
					}
					else {
						print {$client} "CHDIR: ".abs_path("./")."\n";
						print {$client} "NEW: $gene\n";
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

# Don't think this is needed
foreach my $pid (@pids) {
	waitpid($pid, 0);
}

print "\n  All connections closed.\n";
print "Total execution time: ", sec2human(time() - $time), ".\n\n";

# Go back to project directory, delete empty gene dir
chdir("..");
&INT_handler;

sub client {
	my ($opt_name, $address) = @_;	

	my ($server_ip, $port) = split(":", $address);

	chdir("/tmp");
	my $mb = "/tmp/mb";

	#my $pgrp = getpgrp();
	my $pgrp = $$;
	setpgrp();

	# Determine this host's IP
	chomp(my $ip = `dig +short myip.opendns.com \@resolver1.opendns.com`); 

	# Set IP to localhost if we don't have internet
	if ($ip !~ /(?:[0-9]{1,3}\.){3}[0-9]{1,3}/) {
		$ip = "127.0.0.1";
	}

	# Spawn more clients
	my @pids;
	my $total_forks = get_free_cpus(); 
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

	# The name of the gene we are working on
	my $gene;

	# Stores filenames of unneeded files
	my @unlink;

	# Change signal handling so killing the server kills these processes and cleans up
	$SIG{CHLD} = 'IGNORE';
	$SIG{HUP}  = sub { unlink($0, $mb); kill -15, $$; };
	$SIG{TERM} = sub { unlink(glob($gene."*")) if defined($gene); exit(0)};

	# Connect to the server
	my $sock = new IO::Socket::INET(
		PeerAddr  => $server_ip.":".$port,
		Proto => 'tcp') 
	or exit(0); 
	$sock->autoflush(1);

	print {$sock} "NEW: $ip\n";
	while (chomp(my $response = <$sock>)) {

		if ($response =~ /SEND_FILE: (.*)/) {
			my $file_name = $1;
			receive_file({'FILE_PATH' => $file_name, 'FILE_HANDLE' => $sock});	
		}
		elsif ($response =~ /CHDIR: (.*)/) {
			chdir($1);
		}
		elsif ($response =~ /NEW: (.*)/) {
			$gene = $1;

			# Redirect STDOUT to a log file
			open(my $std_out, ">&", *STDOUT);
			open(STDOUT, ">", $gene.".log");

			system($mb, $gene);

			# Put STDOUT back to normal
			open(STDOUT, ">&", $std_out);
			close($std_out);

			unlink($gene);

			# Zip and tarball the results
			my @results = glob($gene."*");
			my $gene_archive_name = "$gene.tar.gz";
			@results = grep {!/\Q$gene_archive_name\E/} @results;
			system("tar", "czf", $gene_archive_name, @results);
			unlink(@results);

			# Send the results back to the server if this is a remote client
			if ($server_ip ne "127.0.0.1" && $server_ip ne $ip) {
				send_file({'FILE_PATH' => $gene_archive_name, 'FILE_HANDLE' => $sock});	
				unlink($gene_archive_name);
			}

			# Request a new job
			print {$sock} "DONE $gene_archive_name || NEW: $ip\n";
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
		unlink($0, $mb);
	}

	exit(0);
}

sub check_nonconvergent {
	my ($opt_name, $threshold) = @_;	

	# We have to do weird things here to get the input name

	my @ARGV = split(/\s+/, $invocation);
	shift(@ARGV); shift(@ARGV);

	# Look for a directory in arguments provided
	my $archive;
	foreach my $arg (@ARGV) {
		if (-d $arg) {
			$archive = $arg;	
		}
	}

	# Die if user didn't give us a directory
	if (!defined($archive)) {
		print "You must specify a directory previously generated by this script to check for nonconvergent genes.\n";
		exit(0);
	}

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

	chdir($project_name);

	# Extract name information from input file
	(my $archive_root = $archive) =~ s/.*\/(.*)/$1/;
	(my $archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar\.gz)|(\.tgz)/$2/;

	# Should have some completed genes in it
	my $incomplete_archive = $archive_root_no_ext.".mb.tar";

	# Check that the incomplete archive exists
	if (!-e $incomplete_archive) {
		print "Could not locate an archive containing completed MrBayes runs.\n";
		exit(0);
	}

	print "\nScript was called as follows:\n$invocation\n\n";

	# Create a temporary directory for our operations
	my $check_dir = "tmp/";
	mkdir($check_dir) if (!-e $check_dir);

	$SIG{INT} = sub { remove_tree($check_dir); exit(0) };

	# Open tarball in genes directory
	chomp(my @genes = `tar xvf '$incomplete_archive' -C $check_dir 2>&1`);
	@genes = map { s/x //; $_ } @genes if ($os_name eq "darwin");
	@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
					(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
					$a <=> $b } @genes;
	my $longest_name_length = length($genes[$#genes]);
	
	print "MrBayes results available for ", scalar(@genes), " total genes:\n";

	chdir($check_dir);

	$SIG{INT} = sub { chdir(".."); remove_tree($check_dir); exit(0) };

	# Parse log of each gene to determine final standard deviation of split frequencies

	my $count = 0;
	foreach my $gene (@genes) {

		chomp(my @contents = `tar xvf '$gene' 2>&1`);
		@contents = map { s/x //; $_ } @contents if ($os_name eq "darwin");

		(my $log_file_path = $gene) =~ s/\.tar\.gz$/.log/;

		# Check log file exists
		if (!-e $log_file_path) {
			print "Could not locate log file for '$gene'.\n";
			exit(0);
		}

		open(my $log_file, "<", $log_file_path);
		chomp(my @data = <$log_file>);
		close($log_file);

		my @splits = grep { /Average standard deviation of split frequencies:/ } @data;
		my $final_split = pop(@splits);

		$final_split =~ s/.*frequencies: (.*)/$1/;

		#print "  $gene: $final_split\n";
		printf("  %-${longest_name_length}s: %s\n", $gene, $final_split);

		if (!defined($final_split) || $final_split > $threshold) {
			$count++;
		}
		unlink(@contents);
	}
	printf("%d gene(s) failed to meet the threshold of %s (%.2f%%).\n", $count, $threshold, ($count / scalar(@genes) * 100));

	# Clean up and exit
	kill(2, $$);
}

sub remove_nonconvergent {
	my ($opt_name, $threshold) = @_;	

	# We have to do weird things here to get the input name

	my @ARGV = split(/\s+/, $invocation);
	shift(@ARGV); shift(@ARGV);

	# Look for a directory in arguments provided
	my $archive;
	foreach my $arg (@ARGV) {
		if (-d $arg) {
			$archive = $arg;	
		}
	}

	# Die if user didn't give us a directory
	if (!defined($archive)) {
		print "You must specify a directory previously generated by this script to check for nonconvergent genes.\n";
		exit(0);
	}

	my $initial_archive = $archive;

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

	chdir($project_name);

	# Extract name information from input file
	(my $archive_root = $archive) =~ s/.*\/(.*)/$1/;
	(my $archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar\.gz)|(\.tgz)/$2/;

	# Should have some completed genes in it
	my $incomplete_archive = $archive_root_no_ext.".mb.tar";

	# Check that the incomplete archive exists
	if (!-e $incomplete_archive) {
		print "Could not locate an archive containing completed MrBayes runs.\n";
		exit(0);
	}

	print "\nScript was called as follows:\n$invocation\n\n";

	# Create a temporary directory for our operations
	my $check_dir = "tmp/";
	mkdir($check_dir) if (!-e $check_dir);

	$SIG{INT} = sub { remove_tree($check_dir); exit(0) };

	# Open tarball in genes directory
	chomp(my @genes = `tar xvf '$incomplete_archive' -C $check_dir 2>&1`);
	@genes = map { s/x //; $_ } @genes if ($os_name eq "darwin");
	@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
					(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
					$a <=> $b } @genes;
	my $longest_name_length = length($genes[$#genes]);
	
	print "MrBayes results available for ", scalar(@genes), " total genes:\n";

	chdir($check_dir);

	$SIG{INT} = sub { chdir(".."); remove_tree($check_dir); exit(0) };

	# Parse log of each gene to determine final standard deviation of split frequencies

	my $count = 0;
	foreach my $gene (@genes) {

		chomp(my @contents = `tar xvf '$gene' 2>&1`);
		@contents = map { s/x //; $_ } @contents if ($os_name eq "darwin");

		(my $log_file_path = $gene) =~ s/\.tar\.gz$/.log/;

		# Check log file exists
		if (!-e $log_file_path) {
			print "Could not locate log file for '$gene'.\n";
			exit(0);
		}

		open(my $log_file, "<", $log_file_path);
		chomp(my @data = <$log_file>);
		close($log_file);

		my @splits = grep { /Average standard deviation of split frequencies:/ } @data;
		my $final_split = pop(@splits);

		$final_split =~ s/.*frequencies: (.*)/$1/;

		#print "  $gene: $final_split";
		if (!defined($final_split) || $final_split > $threshold) {
			unlink($gene);
			#print " -- REMOVED\n";
			printf("  %-${longest_name_length}s: %s -- REMOVED\n", $gene, $final_split);
			$count++;
		}
		else {
			printf("  %-${longest_name_length}s: %s\n", $gene, $final_split);
			#print "\n";
		}
		unlink(@contents);
	}
	printf("%d gene(s) failed to meet the threshold of %s (%.2f%%) and have been removed.\n", $count, $threshold, ($count / scalar(@genes) * 100));

	# Determine which genes met threshold and still remain
	@genes = glob($archive_root_no_ext."*.nex.tar.gz");
	@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
					(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
					$a <=> $b } @genes;

	# Recreate archive with remaining genes
	if (@genes) {
		system("tar", "cf", $incomplete_archive, @genes);
		unlink(@genes);
		system("mv", $incomplete_archive, "..");
	}
	else {
		# Delete the working directory if no genes meet the threshold
		print "No genes met the threshold, removing specified directory.\n";

		chdir($initial_directory);
		$SIG{INT} = sub { remove_tree($initial_archive); exit(0) };
	}

	# Clean up and exit
	kill(2, $$);
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
		#kill(-1, $pid);
		kill(15, $pid);
	}

	# Move into gene directory
	#chdir("$initial_directory");
	chdir("$initial_directory/$project_name");

	# Try to delete directory five times, if it can't be deleted print an error message
	# I've found this method is necessary for analyses performed on AFS drives
	my $count = 0;
	until (!-e $gene_dir || $count == 5) {
		$count++;

		remove_tree($gene_dir, {error => \my $err});
		sleep(1);
	}
	#logger("Could not clean all files in './$gene_dir/'.") if ($count == 5);
	print "Could not clean all files in './$gene_dir/'.\n" if ($count == 5);

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

	return $max_forks if (defined($max_forks));

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

sub usage {
	return "Usage: mb.pl ([PARTITION TARBALL] [-m MRBAYES BLOCK]) || ([MRBAYES TARBALL] [-c THRESHOLD] || [-r THRESHOLD])\n";
}

sub help {
print <<EOF; 
@{[usage()]}
Parallel execution of MrBayes on a large dataset

  -m, --mb-block      text file containing MrBayes commands to append to each input partition (REQUIRED)
  -c, --check         outputs how many MrBayes runs standard deviation of split frequencies have reached the specified threshold
  -r, --remove        removes MrBayes runs with standard deviation of split frequencies below the specified threshold
  -o, --out-dir       name of the directory to store output files in (default: "mb-" + Unix time of script invocation)
  -T, --n-threads     the number of forks ALL hosts running analyses can use concurrently (default: current number of free CPUs)
  --machine-file      file name containing hosts to ssh onto and perform analyses on, passwordless login MUST be enabled
                      for each host specified in this file
  --port              specifies the port to utilize on the server (Default: 10002)
  -h, --help          display this help and exit
  --usage             display proper script invocation format

Examples:
  perl mb.pl align.tgz -m bayes.txt --machine-file hosts.txt     runs MrBayes on each partition stored in align.tgz using parameters 
                                                                 stored in bayes.txt on computers specified in hosts.txt
  perl mb.pl align.mb.tar --check 0.02                           prints which genes in align.mb.tar have MCMC chains which reached a 
                                                                 standard deviation of split frequencies below 0.02
  perl mb.pl align.mb.tar --remove 0.05                          removes genes in align.mb.tar that have MCMC chains which did not reach
                                                                 a standard deviation of split frequencies of 0.05

Mail bug reports and suggestions to <noah.stenz.github\@gmail.com>
EOF
exit(0);
}
