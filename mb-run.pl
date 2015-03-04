#!/usr/bin/perl
use strict;
use warnings;
#use Cwd;
use Cwd qw(abs_path);
use POSIX;
use IO::Select;
use IO::Socket;
use Digest::MD5;
use Time::HiRes qw(time);
use Getopt::Long;

# Turn on autoflush
$|++;

# Server port
my $port = 10002;

# Stores executing machine hostnames
my @machines;

# Path to text file containing computers to run on
my $machine_file_path;

# MrBayes block which will be used for each run
my $mb_block;

# Where this script is located 
my $script_path = abs_path($0);

# General script settings
my $no_forks;

# How the script was called
my $invocation = "perl mb-run.pl @ARGV";

# Name of output directory
#my $project_name = "alignment-breakdown-".time();
my $project_name = "mb-run-dir";

# Read commandline settings
GetOptions(
	"no-forks"          => \$no_forks,
	"mrbayes-block|m:s" => \$mb_block,
	"machine-file:s"    => \$machine_file_path,
	"check|c:s"         => \&check_nonconvergent,
	"remove|r:s"        => \&remove_nonconvergent,
	"client:s"          => \&client, # for internal usage only
	"help|h"            => sub { print &help; exit(0); },
	"usage"             => sub { print &usage; exit(0); },
);

# Get paths to required executables
my $mb = check_path_for_exec("mb");

#my $archive_path = shift(@ARGV);
my $archive = shift(@ARGV);

# Some error checking
die "You must specify an archive file.\n\n", &usage if (!defined($archive));
die "Could not locate '$archive', perhaps you made a typo.\n" if (!-e $archive);
die "You specified a MrBayes run archive instead of an MDL gene archive.\n" if ($archive =~ /-mb\.tar\.gz$/);
die "You must specify a file containing a valid MrBayes block which will be appended to each gene.\n\n", &usage if (!defined($mb_block));
die "Could not locate '$mb_block', perhaps you made a typo.\n\n" if (!-e $mb_block);

# Determine which machines we will run the analyses on
if (defined($machine_file_path)) {
	print "Fetching machine names listed in '$machine_file_path'...\n";
	open(my $machine_file, '<', $machine_file_path);
	chomp(@machines = <$machine_file>);
	close($machine_file);
	print "  $_\n" foreach (@machines);
}

print "\nScript was called as follows:\n$invocation\n";

# Load MrBayes block into memory
open(my $mb_block_file, "<", $mb_block) or die "Could not open '$mb_block': $!.\n";
chomp(my @mb_block = <$mb_block_file>);
close($mb_block_file);

$mb_block = join("\n", @mb_block);

# Extract name information from input file
(my $archive_root = $archive) =~ s/.*\/(.*)/$1/;
#(my $archive_root_no_ext = $archive) =~ s/.*\/(.*)\..*/$1/;
(my $archive_root_no_ext = $archive) =~ s/(.*\/)?(.*)(\.tar\.gz)|(\.tgz)/$2/;

# Initialize working directory
# Remove conditional eventually
mkdir($project_name) || die "Could not create '$project_name'$!.\n" if (!-e $project_name);
my $archive_abs_path = abs_path($archive);
# Remove conditional eventually
run_cmd("ln -s $archive_abs_path $project_name/$archive_root") if (! -e "$project_name/$archive_root");

chdir($project_name);

# Change how Ctrl+C is interpreted to allow for clean up
$SIG{'INT'} = 'INT_handler';

# Define and initialize directories
#my $gene_dir = $alignment_root."mdl-genes/";
my $gene_dir = "genes/";
mkdir($gene_dir) or die "Could not create '$gene_dir': $!.\n" if (!-e $gene_dir);

# TODO: figure out how to rework rerunning, probably just specify output directory
my %complete_genes;
#my $mb_archive = "$archive_root_no_ext-mb.tar.gz";
my $mb_archive = "$archive_root_no_ext.mb.tar.gz";
##my $mb_archive = "$alignment_name-mb.tar.gz";
#if (-e $mb_archive) {
#	print "\nArchive containing completed MrBayes runs found for this dataset found in '$mb_archive'.\n";
#	print "Completed runs contained in this archive will be removed from the job queue.\n";
#	chomp(my @complete_genes = `tar tf $mb_archive`);
#	foreach my $gene (@complete_genes) {
#		$gene =~ s/\.tar\.gz//;
#		$complete_genes{$gene}++;
#	}
#}

# Clean up files remaining from previous runs, needed?
clean_up({'DIRS' => 0});

chomp(my @genes = `tar xvf $archive -C $gene_dir`);

chdir($gene_dir);

#@genes = glob($alignment_name."*");
@genes = glob($archive_root_no_ext."*.nex");
@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
				(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
				$a <=> $b } @genes;

if (%complete_genes) {
	foreach my $index (reverse(0 .. $#genes)) {
		if (exists($complete_genes{$genes[$index]})) {
			splice(@genes, $index, 1);
		}
	}
}

die "\nAll jobs have already completed.\n\n" if (!@genes);

# Append given MrBayes block to the end of each gene
print "\nAppending MrBayes block to each gene... ";
foreach my $gene (@genes) {
	open(my $gene_file, ">>", $gene) or die "Could not open '$gene': $!.\n";
	print {$gene_file} "\n$mb_block\n";
	close($gene_file);
}
print "done.\n\n";

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

		# Send MrBayes executable to the machine
		system("scp", "-q", $mb, $machine.":/tmp");

		# Execute this perl script on the given machine
		# -tt forces pseudo-terminal allocation and lets us stop remote processes
		exec("ssh", "-tt", "$machine", "perl", "/tmp/$script_name", "--client=$server_ip");
		exit(0);
	}
	else {
		push(@pids, $pid);
	}
}

#chdir($gene_dir);

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
my $num_digits = get_num_digits({'NUMBER' => scalar(@genes)});

# Begin the server's job distribution
while ((!defined($total_connections) || $closed_connections != $total_connections) || $total_connections < $starting_connections) {
	# Contains handles to clients which have sent information to the server
	my @clients = $select->can_read(0);

	# Handle each ready client individually
	CLIENT: foreach my $client (@clients) {

		# Client requesting new connection
		if ($client == $sock) {
			$total_connections++;
			$select->add($sock->accept());
		}
		else {

			# Get client's message
		#	my $response = <$client>;
		#	next if (not defined($response)); # a response should never actually be undefined
			my $response = $client->getline();
			next if (not defined($response)); # a response should never actually be undefined

			if ($response =~ /SEND_FILE: (.*)/) {
				my $file_name = $1;
				receive_file({'FILE_PATH' => $file_name, 'FILE_HANDLE' => $client});	
			}
			if ($response =~ /DONE \|\|/) {
				$complete_count++;				
				printf("  Analyses complete: %".$num_digits."d/%d.\r", $complete_count, scalar(@genes));
			}
			if ($response =~ /NEW: (.*)/) {
				my $client_ip = $1;
				if ($job_number < scalar(@genes)) {
					printf("\n  Analyses complete: %".$num_digits."d/%d.\r", 0, scalar(@genes)) if ($job_number == 0);

					my $gene = $genes[$job_number];

					# Check whether the client is remote or local, send it needed files if remote
					if ($client_ip ne $server_ip) {

						$SIG{CHLD} = 'IGNORE';

						# Fork to perform the file transfer and prevent stalling the server
						my $pid = fork(); 
						if ($pid == 0) {
							send_file({'FILE_PATH' => $gene, 'FILE_HANDLE' => $client});						
							unlink($gene);

							print {$client} "NEW: $gene\n";
							exit(0);
						}
					}
					else {
						print {$client} "CHDIR: ".abs_path("./")."\n";
						print {$client} "NEW: $gene\n";
					}
					$job_number++;
				}
				else {
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
print "Total execution time: ", secs_to_readable({'TIME' => time() - $time}), "\n\n";

die;

my @tarballs = map { local $_ = $_; $_ .= ".tar.gz"; $_ } @genes;

print "Compressing and archiving results... ";
#if (-e $alignment_root.$mb_archive) {
if (-e $archive_root_no_ext.$mb_archive) {
	#chdir($alignment_root);
	chdir("..");
	#system("tar", "xf", $mb_archive, @tarballs);
	#chomp(my @complete_genes = `tar xvf $mb_archive -C $gene_dir`);
	chomp(my @complete_genes = `tar xvf $mb_archive -C $gene_dir`);
	push(@tarballs, @complete_genes);
	chdir($gene_dir);
}
system("tar", "czf", $mb_archive, @tarballs, "--remove-files");
#system("mv", $mb_archive, $alignment_root);
system("mv", $mb_archive, $archive_root_no_ext);
print "done.\n\n";

sub client {
	my ($opt_name, $server_ip) = @_;	
	#chdir($ENV{'HOME'});

	chdir("/tmp");
	my $mb = "/tmp/mb";

	#my $pgrp = getpgrp();
	my $pgrp = $$;

	chomp(my $ip = `dig +short myip.opendns.com \@resolver1.opendns.com`); 
	die "Could not establish an IP address for host.\n" if (not defined $ip);

	my @pids;
	my $total_forks = get_free_cpus(); 
	$total_forks-- if ($ip eq $server_ip); # the server also uses 100% of a cpu
	#my $total_forks = 5;

	my $fork_id = 1;
	if ($total_forks > 1) {
		foreach my $fork (1 .. $total_forks - 1) {

			my $pid = fork();
			if ($pid == 0) {
				last;
			}
			else {
				push(@pids, $pid);
				$fork_id++;
			}
		}
	}

	$SIG{CHLD} = 'IGNORE';
	$SIG{HUP} = sub { unlink($0, $mb); kill -15, $$; exit(0); };

	my $gene;

	my @unlink;
	#$SIG{INT} = sub { unlink(@unlink) };
	#$SIG{INT} = sub { unlink(glob($gene."*")) if defined($gene); unlink($check_file) };
	#$SIG{INT} = sub { unlink(glob($gene."*")) if defined($gene) };
	#$SIG{INT} = sub { unlink(glob($gene."*"), $check_file) if defined($gene) };
	$SIG{TERM} = sub { unlink(glob($gene."*")) if defined($gene); };

	my $sock = new IO::Socket::INET(
		PeerAddr  => $server_ip.":".$port,
		Proto => 'tcp') 
	or exit(0); 
	$sock->autoflush(1);
	#$sock->autoflush(0);

	print {$sock} "NEW: $ip\n";
	#print {$sock} "FILE: $check_file || NEW: $ip\n";
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

			open(my $std_out, ">&", *STDOUT);
			open(STDOUT, ">", $gene.".log");

			system($mb, $gene);

			open(STDOUT, ">&", $std_out);
			close($std_out);

			unlink($gene);

			my @results = glob($gene."*");

			my $gene_archive_name = "$gene.tar.gz";
			system("tar", "czf", $gene_archive_name, @results);

			if ($server_ip ne $ip) {
				send_file({'FILE_PATH' => $gene_archive_name, 'FILE_HANDLE' => $sock});	
				unlink($gene_archive_name);

#				my $pid = fork();
#				if ($pid == 0) {
#					send_file({'FILE_PATH' => $gene_archive_name, 'FILE_HANDLE' => $sock});	
#					unlink($gene_archive_name);
#					unlink(@results);
#
#					print {$sock} "DONE || NEW: $ip\n";
#					exit(0);
#				}
			}
#			else {
#				unlink(@results);
#			}

			unlink(@results);

			print {$sock} "DONE || NEW: $ip\n";
		}
		elsif ($response eq "HANGUP") {
			last;
		}
		$sock->flush();
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

#sub check_nonconvergent {
#	my ($opt_name, $threshold) = @_;	
#
#	print "\nScript was called as follows:\n$invocation\n";
#
#	die "You must specify an archive file.\n\n", &usage if (!defined($archive));
#	die "Could not locate '$archive', perhaps you made a typo.\n\n" if (!-e $archive);
#
#	($alignment_root = $archive) =~ s/(.*\/).*/$1/;
#	($alignment_name = $archive) =~ s/.*\/(.*)-genes\.tar\.gz/$1/;
#	$alignment_name =~ s/.*\/(.*)-mb\.tar\.gz/$1/;
#
#	if ($alignment_root eq $alignment_name) {
#		$alignment_root = getcwd()."/";
#		$alignment_name =~ s/(.*)-genes\.tar\.gz/$1/;
#		$alignment_name =~ s/(.*)-mb\.tar\.gz/$1/;
#	}
#
#	my $archive_name = $alignment_name."-mb.tar.gz";
#
#	chdir($alignment_root);
#
#	die "Could not locate an archive containing completed MrBayes runs.\n" if (!-e $alignment_name."-mb.tar.gz");
#
#	print "\nScript was called as follows:\n$invocation\n";
#
#	my $tmp_dir = "mb-tmp/";
#	mkdir($tmp_dir) if (!-e $tmp_dir);
#
#	system("tar", "xf", $archive_name, "-C", $tmp_dir);
#
#	chdir($tmp_dir);
#	my @genes = glob($alignment_name."*.nex.tar.gz");
#	@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
#					(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
#					$a <=> $b } @genes;
#	
#	print "MrBayes results available for ", scalar(@genes), " total genes.\n";
#	
#	my $count = 0;
#	foreach my $gene (@genes) {
#		chomp(my @contents = `tar tf $gene`);
#		system("tar", "xf", $gene);
#
#		(my $log_file_path = $gene) =~ s/\.tar\.gz$/.log/;
#
#		open(my $log_file, "<", $log_file_path);
#		chomp(my @data = <$log_file>);
#		close($log_file);
#
#		my @splits = grep { /Average standard deviation of split frequencies:/ } @data;
#		my $final_split = pop(@splits);
#
#		$final_split =~ s/.*frequencies: (.*)/$1/;
#
#		#if ($final_split > $threshold) {
#		if (!defined($final_split) || $final_split > $threshold) {
##			unlink($gene);
#			$count++;
#		}
#		unlink(@contents);
#	}
#	print "$count gene(s) failed to meet the threshold ($threshold).\n";
#
#	@genes = glob($alignment_name."*.nex.tar.gz");
##	@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
##					(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
##					$a <=> $b } @genes;
#	unlink(@genes);
#
##	if (@genes) {
##		system("tar", "czf", $archive_name, @genes);
##		system("mv", $archive_name, $alignment_root);
##		unlink(@genes);
##	}
#
#	chdir($alignment_root);
##	unlink($archive_name) if (!@genes); # delete the archive if no genes meet the threshold
#	
#	rmdir($tmp_dir);
#		
#	exit(0);
#}
#
#sub remove_nonconvergent {
#	my ($opt_name, $threshold) = @_;	
#
#	die "You must specify an archive file.\n\n", &usage if (!defined($archive));
#	die "Could not locate '$archive', perhaps you made a typo.\n\n" if (!-e $archive);
#
#	($alignment_root = $archive) =~ s/(.*\/).*/$1/;
#	($alignment_name = $archive) =~ s/.*\/(.*)-genes\.tar\.gz/$1/;
#	$alignment_name =~ s/.*\/(.*)-mb\.tar\.gz/$1/;
#
#	if ($alignment_root eq $alignment_name) {
#		$alignment_root = getcwd()."/";
#		$alignment_name =~ s/(.*)-genes\.tar\.gz/$1/;
#		$alignment_name =~ s/(.*)-mb\.tar\.gz/$1/;
#	}
#
#	my $archive_name = $alignment_name."-mb.tar.gz";
#
#	chdir($alignment_root);
#
#	die "Could not locate an archive containing completed MrBayes runs.\n" if (!-e $alignment_name."-mb.tar.gz");
#
#	print "\nScript was called as follows:\n$invocation\n";
#
#	my $tmp_dir = "mb-tmp/";
#	mkdir($tmp_dir) if (!-e $tmp_dir);
#
#	system("tar", "xf", $archive_name, "-C", $tmp_dir);
#
#	chdir($tmp_dir);
#	my @genes = glob($alignment_name."*.nex.tar.gz");
#	@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
#					(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
#					$a <=> $b } @genes;
#	
#	print "MrBayes results available for ", scalar(@genes), " total genes.\n";
#	
#	my $count = 0;
#	foreach my $gene (@genes) {
#		chomp(my @contents = `tar tf $gene`);
#		system("tar", "xf", $gene);
#
#		(my $log_file_path = $gene) =~ s/\.tar\.gz$/.log/;
#
#		open(my $log_file, "<", $log_file_path);
#		chomp(my @data = <$log_file>);
#		close($log_file);
#
#		my @splits = grep { /Average standard deviation of split frequencies:/ } @data;
#		my $final_split = pop(@splits);
#
#		$final_split =~ s/.*frequencies: (.*)/$1/;
#		
#		#if ($final_split > $threshold) {
#		if (!defined($final_split) || $final_split > $threshold) {
#			unlink($gene);
#			$count++;
#		}
#		unlink(@contents);
#	}
#	print "$count gene(s) failed to meet the threshold ($threshold) and have been removed.\n";
#
#	@genes = glob($alignment_name."*.nex.tar.gz");
#	@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
#					(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
#					$a <=> $b } @genes;
#
#	if (@genes) {
#		system("tar", "czf", $archive_name, @genes);
#		system("mv", $archive_name, $alignment_root);
#		unlink(@genes);
#	}
#
#	chdir($alignment_root);
#	unlink($archive_name) if (!@genes); # delete the archive if no genes meet the threshold
#	
#	rmdir($tmp_dir);
#		
#	exit(0);
#}

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
		kill(9, $pid);
	}

	chdir($gene_dir);
#	my @genes = glob($alignment_name."*.nex.tar.gz");
#	@genes = sort { (local $a = $a) =~ s/.*-(\d+)-\d+\..*/$1/; 
#					(local $b = $b) =~ s/.*-(\d+)-\d+\..*/$1/; 
#					$a <=> $b } @genes;
#	chdir($alignment_root);
#
#	if (@genes) {
#		print "Adding results from ".scalar(@genes)." run(s) to '$mb_archive'.\n";
#		if (-e $mb_archive) {
#			chomp(my @complete_genes = `tar tf $mb_archive`);
#			system("tar", "xf", $mb_archive, "-C", $gene_dir);
#			push(@genes, @complete_genes);
#			unlink($mb_archive);
#		}
#		chdir($gene_dir);
#		system("tar", "czf", $mb_archive, @genes);
#		system("mv", $mb_archive, $alignment_root);
#	}
#
#	clean_up({'DIRS' => 1});
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

sub help {

}

sub usage {

}
