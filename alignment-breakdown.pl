#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);
use Fcntl;
use Fcntl ':flock';
use POSIX;
use IO::Pipe;
use IO::Handle;
use IO::Select;
use IO::Socket;
use Digest::MD5;
use Getopt::Long;
use Time::HiRes 'time';

# prevents printing a recursion warning which appear with large datasets
no warnings 'recursion';

# Turn on autoflush
$|++;

# Pseudo-boolean that stores IP of server
my $client;

# Server port
my $port = 10001;

# Stores executing machine hostnames
my @machines;

# Path to text file containing computers to run on
my $machine_file_path;

# Where this script is located 
my $script_path = abs_path($0);

# MDL settings
my $gap_isnt_char;
my $exclude_gaps;
my $nletters = 4;
my $nbestpart = 1;
my $ngroupmax = 10000;

my $no_forks;
my $forced_break;
my $min_block_size;

# Stores total number of parsimony-informative characters found by script
my $num_pars_inf_chars;

# How the script was called
my $invocation = "perl alignment-breakdown.pl @ARGV";

# Name of output directory
#my $project_name = "alignment-breakdown-".time();
my $project_name = "alignment-breakdown-dir";

# Get paths to required executables
my $mdl = check_path_for_exec("mdl");
my $paup = check_path_for_exec("paup");

# Read commandline settings
GetOptions(
	"block-size|b:i"   => \$min_block_size,
	"nletters|l:i"     => \$nletters,
	"nbestpart|p:i"    => \$nbestpart,
	"ngroupmax|m:i"    => \$ngroupmax,
	#"gap-as-char|g"    => \$gap_as_char,
	"gap-isnt-char|g"  => \$gap_isnt_char,
	"exclude-gaps|e"   => \$exclude_gaps,
	"forced-break|f:i" => \$forced_break,
	"no-forks"         => \$no_forks,
	"machine-file:s"   => \$machine_file_path,
	"client:s"         => \&client, # for internal usage only
	"help|h"           => sub { print &help; exit(0); },
	"usage"            => sub { print &usage; exit(0); },
);

my $align = shift(@ARGV);

# Some error checking
die "You must specify an alignment file.\n\n", &usage if (!defined($align));
die "Could not locate '$align', perhaps you made a typo.\n" if (!-e $align);
die "You must specify an alignment to partition.\n\n", &usage if (!defined($align));
die "You must specify a minimum block length.\n\n",    &usage if (!defined($min_block_size));
die "The minimum block size must be greater than 0!\n" if ($min_block_size <= 0);
die "Forced breakpoints cannot be placed at the distance you specified!\n" if (defined($forced_break) && $forced_break <= 0);
die "The minimum block size can not be greater than or equal to the forced breakpoint length.\n" if (defined($forced_break) && $min_block_size >= $forced_break);

print "WARNING: It is not recommended to have forced breakpoints placed so frequently unless you are analyzing a very large dataset.\n" if (defined($forced_break) && $forced_break <= 100);

# Adjust MDL settings depending on user input
$gap_isnt_char = 1 if (defined($exclude_gaps));
$nletters++        if (!defined($exclude_gaps));

# Extract name information from input file
(my $align_root = $align) =~ s/.*\/(.*)/$1/;
(my $align_root_no_ext = $align) =~ s/.*\/(.*)\..*/$1/;

# Initialize working directory
mkdir($project_name) || die "Could not create '$project_name'$!.\n" if (!-e $project_name);

my $align_abs_path = abs_path($align);
# Remove conditional eventually
run_cmd("ln -s $align_abs_path $project_name/$align_root") if (! -e "$project_name/$align_root");

# Determine which machines we will run the analyses on
if (defined($machine_file_path)) {
	print "Fetching machine names listed in '$machine_file_path'...\n";
	open(my $machine_file, '<', $machine_file_path);
	chomp(@machines = <$machine_file>);
	close($machine_file);
	print "  $_\n" foreach (@machines);
}

chdir($project_name);

# Define and initialize directories
my $gene_dir      = "mdl-genes/";
my $score_dir     = "mdl-scores/";
my $phylip_dir    = "mdl-phylip/";
my $partition_dir = "mdl-partitions/";

mkdir($gene_dir)      or die "Could not create '$gene_dir': $!.\n"      if (!-e $gene_dir);
mkdir($score_dir)     or die "Could not create '$score_dir': $!.\n"     if (!-e $score_dir);
mkdir($phylip_dir)    or die "Could not create '$phylip_dir': $!.\n"    if (!-e $phylip_dir);
mkdir($partition_dir) or die "Could not create '$partition_dir': $!.\n" if (!-e $partition_dir);

# Change how Ctrl+C is interpreted to allow for clean up
$SIG{'INT'} = 'INT_handler';

# Print the current script settings
print "\nScript was called as follows:\n$invocation\n";

# Output some more useful information on current settings
if (defined($forced_break)) {
	print "\nWill now proceed to breakdown '$align' using a forced breakpoint after ";
	print "",(($forced_break == 1) ? "every character, " : "every $forced_break characters, "), "and a minimum block size of $min_block_size.\n";
}
else {
	print "\nWill now proceed to breakdown '$align' using a minimum block size of $min_block_size.\n";
}
print "MDL settings: nletters = $nletters, nbestpart = $nbestpart, ngroupmax = $ngroupmax.\n\n";

# Relative site numbers of PI characters
my @locations;

# Absolute site numbers of PI characters
my @translated_locations;

#parse_input();
my %align = parse_input();
get_informative_chars({'ALIGN' => \%align});
run_mdl();

sub parse_input {
	my $pipe = new IO::Pipe;
	my $select = new IO::Select; 
	my $pid = fork();

	my %align;
	# By handling the file reading within a fork we save some RAM
	if ($pid == 0) {
		my $TO_PARENT = $pipe->writer();

		# Open the input file
		open(my $align_file, '<', $align_root) 
			or die "Could not open '$align_root': $!\n";
		chomp(my @data = <$align_file>);
		close($align_file);

		# Determine whether file is in Nexus or FASTA format, then parse it
		if ($data[0] =~ /#NEXUS/i) {
			#print "Input file \"$alignment_path\" appears to be a  Nexus file.\n\n";
			print "Input file '$align' appears to be a  Nexus file.\n";
			parse_nexus({'ALIGN' => \%align, 'DATA' => \@data});
		}
		elsif ($data[0] =~ /^>/) {
			#print "Input file \"$alignment_path\" appears to be a FASTA file.\n\n";
			print "Input file '$align' appears to be a FASTA file.\n";
			parse_fasta({'ALIGN' => \%align, 'DATA' => \@data});
		}

		# Send sequence info back to parent
		foreach my $key (keys %align) {
			print {$TO_PARENT} "$key||$align{$key}\n";
		}

		exit(0);
	}
	else {
		my $FROM_CHILD = $pipe->reader();
		$select->add($FROM_CHILD);
	}

	# Reap children
	while (my @handles = $select->can_read()) {
		my @handles = $select->can_read();
		my $handle = shift(@handles);

		chomp(my @child_response = <$handle>);
		foreach my $response (@child_response) {
			my ($key, $value) = split(/\|\|/, $response);
			$align{$key} = $value;
		}
		waitpid($pid, 0);
		$select->remove($handle);
		$handle->close;
	}

	return %align;
}

# Parses alignment data from a properly formatted Nexus file
# Parameters:
#     DATA  - an array reference containing the entire contents of the Nexus file
#     ALIGN - a reference to a hash which will have the alignment data stored with 
#             taxa names as itskeys and the nucleotide data as values
sub parse_nexus {
	my $settings = shift;

	my $align = $settings->{'ALIGN'};

	# TODO: Rewrite this, pretty sure it isn't very good
	
	my $in_data_block = 0;
	my $in_data_matrix = 0;
	foreach my $line (@{$settings->{'DATA'}}) {
		if ($line =~ /Begin data;/i) {
			$in_data_block = 1;
		}
		elsif ($in_data_block && $line =~ /end;/i) {
			$in_data_block = 0;
		}
		elsif ($in_data_block && $line =~ /Matrix/i) {
			$in_data_matrix = 1;
		}
		elsif ($in_data_matrix && $line =~ /;/) {
			$in_data_matrix = 0;
		}
		elsif ($in_data_block && $in_data_matrix) {
			if ($line =~ /(\S+)\s+(\S+)/) {
				my $taxon = $1;
				my $alignment = $2;
				$align->{$taxon} .= $alignment;
			}
		}
	}
}

# Parses alignment data from a properly formatted FASTA file
# Parameters:
#     DATA  - an array reference containing the entire contents of the FASTA file
#     ALIGN - a reference to a hash which will have the alignment data stored with 
#             taxa names as its keys and the nucleotide data as values
sub parse_fasta {
	my $settings = shift;

	my $align = $settings->{'ALIGN'};
	
	my $taxon;
	foreach my $line (@{$settings->{'DATA'}}) {
		if ($line =~ /^>(.*)/) {
			$taxon = $1;
		}
		else {
			$taxon =~ s/-/_/g;
			$align->{$taxon} .= $line;
		}
	}
}

# Writes a single Nexus file consisting of only parsimony-informative characters
# Parameters:
#     TOTAL - the total number of reduced files which will be output
#     ALIGN - a reference to a hash which has the alignment data stored with taxa
#             names as its keys and the nucleotide data as values
#     NUM   - the number of the specific output file which should be written
#     TOKEN - an identification token which will be placed in the output files
sub write_nexus_file_reduced {
	my $settings = shift;

	my $number = $settings->{'NUM'};
	my $align = $settings->{'ALIGN'};
	my $token = $settings->{'TOKEN'};

	# Adjust the partition length of the output segment to prevent overshoot.
	my $forced_break = $forced_break;
	if ($forced_break * ($number + 1) > scalar(@locations)) {
		$forced_break = -1 * (($forced_break * ($number + 1)) - scalar(@locations) - $forced_break);
	}

	# Use sysopen and the O_NDELAY flag for faster output.
	my $file_name = $align_root_no_ext."-$token-".($number+1).".nex";
	unlink($file_name) if (-e $file_name);

	sysopen(my $out, $file_name, O_WRONLY|O_NDELAY|O_CREAT)
		or die "Could not open '$file_name': $!.\n";

	my @taxa = sort { $a cmp $b } (keys %{$align});

	print {$out} "#NEXUS\nbegin data;\n dimensions ntax=".scalar(@taxa);
	print {$out} " nchar=$forced_break;\n format datatype=dna gap=- missing=?;\n matrix\n";

	foreach my $taxon (@taxa) {
		print {$out} " $taxon\t";
		for (my $i = 0; $i < $forced_break; $i++ ) {
			print {$out} substr($align->{$taxon}, $locations[$i + ($number * $forced_break)] - 1, 1);
		}
		print {$out} "\n";
	}	

	# This would be needed in order to allow continuation of a failed run
	# Add the absolute positions as a comment on the bottom for later scripts.

#	print {$out} " ;\nend;\n[Absolute Character Locations:]\n[";
#	for (my $i = 0; $i < $forced_break; $i++) {
#		#print {$out} @{$locations}[$i +  ($number * $forced_break)],"\t";
#		print {$out} $locations[$i +  ($number * $forced_break)],"\t";
#	}
#	print {$out} "]\n";

	print {$out} " ;\nend;\n";
	close($out);
}

# Writes a single Nexus file
# Parameters:
#     TOTAL - the total number of reduced files which will be output
#     ALIGN - a reference to a hash which has the alignment data stored with taxa
#             names as its keys and the nucleotide data as values
#     NUM   - the number of the specific output file which should be written
#     TOKEN - an identification token which will be placed in the output files
sub write_nexus_file {
	my $settings = shift;

	my $number = $settings->{'NUM'};
	my $align = $settings->{'ALIGN'};
	my $total = $settings->{'TOTAL'};
	my $token = $settings->{'TOKEN'};

	# Determine which proportion to write
	my $beginning_fraction = $number / $total;
	my $end_fraction = ($number + 1) / $total;

	my $chromosome_length = length((values %{$align})[0]);
	
	# Determine the actual start and end site numbers
	my $start = int($chromosome_length * $beginning_fraction);
	my $end = int($chromosome_length * $end_fraction);

	my $nchar = $end - $start;

	my $file_name = $align_root_no_ext."-$token-".($number+1)."of$total.nex";
	unlink($file_name) if (-e $file_name);

	sysopen(my $out, $file_name, O_WRONLY|O_NDELAY|O_CREAT)
		or die "Could not open '$file_name': $!.\n";

	my $line_length = 90000;
	my @taxa = sort { $a cmp $b } (keys %{$align});

	print {$out} "#NEXUS\nbegin data;\n dimensions ntax=".scalar(@taxa)." nchar=$nchar;\n ";
	print {$out} "format datatype=dna interleave=yes gap=- missing=?;\n matrix\n";

	for (my $i = 0; $i < $nchar; $i += $line_length) {
		for (my $j = 0; $j < scalar(@taxa); $j++) {
			my $taxon = $taxa[$j];
			my $current_char = $start + $i;
			my $interval = $line_length;

			# Prevent spillover to next block			
			if ($end < $current_char + $interval) {
				$interval = $end - $current_char;
			}
			print {$out} "  ".$taxon."\t".substr($align->{$taxon}, $start + $i, $interval),"\n";
		}
	}
	print {$out} "\n ;\nend;";
	close($out);
}

# Determines which sites are parsimony-informative using PAUP's 'cstatus' command
# Parameters:
#     TOTAL     - the total number of reduced files which will be output
#     NUM       - the number of the specific output file which should be written
#     TOKEN     - an identification token which will be placed in the output files
#     TO_PARENT - a handle used to send information to the parent of the process
sub paup_cstatus {
	my $settings = shift;

	my $number = $settings->{'NUM'};
	my $total = $settings->{'TOTAL'};
	my $token = $settings->{'TOKEN'};
	my $TO_PARENT = $settings->{'TO_PARENT'};

	my $file_name = $align_root_no_ext."-$token-".($number+1)."of$total.nex";

	my $paup_commands = "execute $file_name;\ncstatus full=yes;\n";

	# File name of log which contains PI character status
	my $log_path = "log-$align_root_no_ext"."-$token-".($number+1)."of$total.nex";

	# Temporarily redirect STDOUT to prevent clutter
	open(my $std_out, ">&", *STDOUT);
	close(STDOUT);

	# Run PAUP via a pipe to avoid having to create a command file
	open(my $paup_pipe, "|-", $paup, "-n", "-l", $log_path); 
	foreach my $command (split("\n", $paup_commands)) {
		print {$paup_pipe} $command,"\n";
	}
	close($paup_pipe);

	# Redirect STDOUT back to terminal
	open(STDOUT, ">&", $std_out);
	close($std_out);

	# Read the resulting log with sysread (faster) 
	my $log_contents;
	open(my $log, "<", $log_path) 
		or die "Could not open '$log_path': $!.\n";
	sysread($log, $log_contents, -s $log_path);
	close($log) and unlink($log_path);

	my @log_lines = split("\n", $log_contents);

	# Use map to create an array of only parsimony-informative sites
	my @filtered_log = grep(/^\d+\s+|(Data matrix)/, @log_lines);
	my @locations = map { my @line = split(" ", $_);  ($line[2] eq '-' && $line[3] ne 'X') ? $line[0] : () } @filtered_log;

	# Determine total number of characters in the alignment
	my $total_chars;
	my $info_line = shift(@filtered_log);
	if ($info_line =~ /data matrix has \d+ taxa, (\d+) characters/i) { 
		$total_chars = $1;
	}

	# Send the locations of the parsimony-informative characters to the parent.
	foreach my $location (@locations) {
		print {$$TO_PARENT} "$number-$total_chars--$location\n";
	}

	# Cleanup
	unlink($log_path);
	unlink($file_name);
}

# Determines which sites are parsimony-informative
# Parameters:
#     ALIGN - a reference to a hash which has the alignment data stored with taxa
#             names as its keys and the nucleotide data as values
sub get_informative_chars {
	my $settings = shift;

	my %align = %{$settings->{'ALIGN'}};

	# Split the original input alignment and write a new Nexus file for each split
	my $free_cpus = get_free_cpus();
	parallelize({'METHOD'      => 'write_nexus_file', 
	             'METHOD_ARGS' => {'ALIGN' => \%align, 
				                   'TOKEN' => 'unreduced'}, 
				 'MAX_FORKS'   => $free_cpus, 
				 'TOTAL'       => $free_cpus});

	# Run the PAUP cstatus command to determine which characters are parsimony-informative
	my @raw_locations = parallelize_with_return({'METHOD'      => 'paup_cstatus', 
	                                             'METHOD_ARGS' => {'TOKEN' => 'unreduced'}, 
												 'MAX_FORKS'   => $free_cpus, 
												 'TOTAL'       => $free_cpus});

	# Preparation for determining PIC site numbers
	my %locations;
	foreach my $response (@raw_locations) {
		if ($response =~ /(\d+-\d+)--(\d+)$/) {
			push(@{$locations{$1}}, $2);
		}
	}

	# Calculate absolute PIC sites
	my $start = 0;
	foreach my $key (sort { (local $a = $a) =~ s/(\d+).*/$1/; (local $b = $b) =~ s/(\d+).*/$1/; $a <=> $b } keys %locations) {
		foreach my $position (@{$locations{$key}}) {
			push(@locations, $position + $start);	
		}
		(my $characters = $key) =~ s/(\d+)-//;
		$characters =~ s/--.*//;
		$start += $characters;
	}
	undef(%locations);

	$num_pars_inf_chars = scalar(@locations);

	print "\n", $num_pars_inf_chars, " total parsimony-informative sites found for '$align'.\n\n";

	# Determine how many times the alignment will be forcibly segmented
	my $total_output_files = 1;
	if (defined($forced_break)) {
		$total_output_files = ceil($num_pars_inf_chars / $forced_break);
	}
	else {
		$forced_break = $num_pars_inf_chars;
	}

	# Output the reduced Nexus files containing only PICs
	parallelize({'METHOD' => 'write_nexus_file_reduced', 
	             'METHOD_ARGS' => {'ALIGN' => \%align, 
				                   'TOKEN' => 'reduced', 
								   'LOCATIONS' => \@locations}, 
                 'MAX_FORKS' => $free_cpus, 
				 'TOTAL' => $total_output_files});
}

# Runs the actual MDL analyses required to breakdown the given alignment
sub run_mdl {
	my @genes;
	my %break_info;

	# Determine how many files the reduced sites were split into
	my $total_output_files = 1;
	if (defined($forced_break)) {
		$total_output_files = ceil($num_pars_inf_chars / $forced_break);
	}
	else {
		$forced_break = $num_pars_inf_chars;
	}

	# Perform some preprocessing on the partitions to make later calculations easier

	my %nchar;
	my %genes;
	my $count = 1;
	my $total = ceil($num_pars_inf_chars / $forced_break);
	for (my $i = 0; $i < $num_pars_inf_chars; $i += $forced_break) {

		my $align_length = $forced_break;
		if ($i + $forced_break > $num_pars_inf_chars) {
			$align_length = $num_pars_inf_chars - (($count - 1) * $forced_break);
		}
		$nchar{$count} = $align_length;
		
		my $input_file_name = $align_root_no_ext."-reduced-".$count."of"."$total.nex";

		my $last_gene_size = $align_length % $min_block_size;
		my $total_genes = ($align_length - $last_gene_size) / $min_block_size;
		$total_genes++ if ($last_gene_size);

		my $total_blocks = get_num_blocks({'TOTAL_GENES'  => $total_genes, 
										   'FILE_NUMBER'  => $count, 
		                                   'ALIGN_LENGTH' => $align_length});

		$break_info{$count}->{'TOTAL_GENES'} = $total_genes;
		$break_info{$count}->{'ALIGN_LENGTH'} = $align_length;
		$break_info{$count}->{'TOTAL_BLOCKS'} = $total_blocks;

		$count++;
	}

	chdir($project_name);

	# Remove any files that may remain from previous runs
	clean_up({'DIRS' => 0});

	# Determine how many PAUP iterations must be performed
	my $total_blocks = 0; 
	foreach my $key (keys %break_info) {
		$total_blocks += $break_info{$key}->{'TOTAL_BLOCKS'};
	}

	print "Parsimony analyses will be performed on $total_blocks different blocks.\n";

	# Print a warning if this analysis is going to take a while
	if ($total_blocks > 1000000) {
		print "This analysis will require running more than 1 million parsimony calculations to finish. Okay to continue? [y/n]: ";
		while (chomp(my $input = <STDIN>)) {
			last if ($input =~ /^y/i);
			exit if ($input =~ /^n/i);
			print "[y/n]: ";
		}
		close(STDIN);
	}
	print "\n";

	# TODO: more reliable way to get external IP
	#chomp(my $server_ip = `wget -qO- ipecho.net/plain`); # Returns the external IP address of this computer
	chomp(my $server_ip = `dig +short myip.opendns.com \@resolver1.opendns.com`); # Returns the external IP address of this computer
	die "Could not establish the IP address for the server.\n" if (!defined $server_ip);

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

	# Determine server hostname and add to machines if none were specified by the user
	chomp(my $server_hostname = `hostname`);
	push(@machines, $server_hostname) if (scalar(@machines) == 0);

	my @pids;
	foreach my $machine (@machines) {

		# Fork and create a client on the given machine
		my $pid = fork();	
		if ($pid == 0) {
			(my $script_name = $script_path) =~ s/.*\///;

			# Send this script to the machine
			system("scp", "-q", $script_path, $machine.":/tmp");

			# Send PAUP executable to the machine
			system("scp", "-q", $paup , $machine.":/tmp");

			# Send reduced alignments to remote machines
			if ($machine ne $server_hostname) {
				system("scp -q *-reduced-* $machine.:/tmp");
			}

			# Execute this perl script on the given machine
			exec("ssh", $machine, "perl", "/tmp/".$script_name, "--client=$server_ip");

			exit(0);
		}
		else {
			push(@pids, $pid);
		}
	}

	my $select = IO::Select->new($sock);

	# Stores which job is next in queue 
	my $job_number = 0;

	# Number of open connections to a client
	my $total_connections;

	# Number of connections server has closed
	my $closed_connections = 0;

	# Minimum number of connections server should expect
	my $starting_connections = scalar(@machines);

	# Stores parsimiony score iformation for blocks
	my %scores;
	my $time = time();

	# Jobs are delegated in $subset_size queues
	my @gene_subset;
	my $subset_count = 0;
	my $subset_size = 3000;

	# Begin the server's job distribution
	while ((!defined($total_connections) || $closed_connections != $total_connections) || $total_connections < $starting_connections) {

		# Contains handles to clients which have sent information to the server
		my @clients = $select->can_read(0);

		# To prevent %scores from becoming too large, we dump it after every $subset_size calculations
		my $total_scores = 0;
		foreach my $key (keys %scores) {
			$total_scores += scalar(@{$scores{$key}});

			# We have enough scores stored to warrant output
			if ($total_scores >= $subset_size) {
				foreach my $index (sort { $a <=> $b } keys %scores) {

					my @partition = @{$scores{$index}};
					open(my $joined_score_file, '>>', $score_dir.$align_root_no_ext."-all-scores-$index") 
						or die "Could not open '".$score_dir.$align_root_no_ext."-all-scores-$index': $!.\n";
					foreach my $line (@partition) {
						print {$joined_score_file} $line,"\n";	
					}	
					close($joined_score_file);
				}
				undef(%scores);
			}
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

				# Client has finished a job
				if ($response =~ /RESULT: (.*):(\S+)/) {

					my ($file_name, $score) = ($1, $2);
					(my $partition = $file_name) =~ s/.*-scores-(\d+)-\d+.*/$1/;
					(my $tree_number = $file_name) =~ s/.*-scores-\d+-(\d+).*/$1/;

					push(@{$scores{$partition}}, "$tree_number : $score");
					if ($job_number % 100 == 0 || $job_number == $total_blocks) {
						my $num_digits = get_num_digits({'NUMBER' => $total_blocks});
						printf("    Analyses complete: %".$num_digits."d/%d.\r", $job_number, $total_blocks);
					}

				}

				# Client wants a new job
				if ($response =~ /NEW: (.*)/) {

					my $client_ip = $1;

					# Check that we have jobs to delegate
					if ($job_number < $total_blocks) {

						# Check if our current queue has jobs in it
						if (scalar(@gene_subset) == 0) {

							# Determine which forced partition the current job in the queue is in
							my $sum;
							my $break_number;
							foreach my $key (sort { $a <=> $b } keys %break_info) {
								$sum += $break_info{$key}->{'TOTAL_BLOCKS'};
								if ($sum >= $job_number) {
									$break_number = $key;
									last;
								}
							}

							my $start = $job_number;

							# Get a new gene subset
							@gene_subset = get_block_subset({'FILE_NUMBER'  => $break_number, 
															 'START' => $start,
															 'SIZE' => $subset_size,
															 'GENES' => \@gene_subset,
															 'BREAK_INFO' => \%break_info});

							# Write the command files required to run each job
							print "\n  Writing command files for each block... ";
							
							# TODO: eliminate actual writing of command, have it return file contents
							my $time = time();
							my $free_cpus = get_free_cpus();
					#		parallelize({'METHOD'      => 'write_mdl_command', 
					#					 'METHOD_ARGS' => {'GENES' => \@gene_subset}, 
					#					 'MAX_FORKS'   => $free_cpus, 
					#					 'TOTAL'       => scalar(@gene_subset)});
							write_mdl_command({'GENES' => \@gene_subset});

							print "done. ($start - ".($start + scalar(@gene_subset) - 1).")\n";
							print "  Total execution time: ", secs_to_readable({'TIME' => time() - $time}), "\n";

							$subset_count++;
						}

						
						# Retrieve the current job and information from the queue
						my $job = shift(@gene_subset);
						my $paup_command = $job->{'PAUP_COMMAND'};
						my $score_file_name = $job->{'SCORE_FILE_NAME'};
						my $align_file_name = $job->{'ALIGN_FILE_NAME'};

						# Check whether the client is remote or local, send it needed files if remote
						if ($client_ip ne $server_ip) {

							$SIG{CHLD} = 'IGNORE';

							# Fork to perform the file transfer and prevent stalling the server
							my $pid = fork(); 
							if ($pid == 0) {
								print {$client} "NEW: '$score_file_name' '$paup_command'\n";
								exit(0);
							}
						}
						else {
							print {$client} "CHDIR: ".abs_path($score_dir)."\n";
							print {$client} "NEW: '$score_file_name' '$paup_command'\n";
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

	# Output any block scores that haven't been dumped yet
	if (scalar(keys %scores) > 0) {
		foreach my $index (sort { $a <=> $b } keys %scores) {
			my @partition = @{$scores{$index}};
			open(my $joined_score_file, '>>', $score_dir.$align_root_no_ext."-all-scores-$index") 
				or die "Could not open '".$score_dir.$align_root_no_ext."-all-scores-$index': $!.\n";
			foreach my $line (@partition) {
				print {$joined_score_file} $line,"\n";	
			}	
			close($joined_score_file);
		}
		undef(%scores);
	}

	print "\n  All connections closed.\n\n";
	print "Total execution time: ", secs_to_readable({'TIME' => time() - $time}), "\n\n";

	chdir($project_name);
	#parse_input({'GET_ALIGN' => 1});
	#get_informative_chars({'ALIGN' => \%align});

	my $ntax = scalar(keys %align);

	# Run each partition defined by the forced break points (if any) through MDL
	foreach my $partition (1 .. $total) {
		my $data_file = "$align_root_no_ext-reduced-$partition"."of$total.nex";
		my $output_name = $partition_dir."$align_root_no_ext-$partition"."of$total.partitions";
		system("$mdl -ntax $ntax -nchar $nchar{$partition} -scorefile mdl-scores/$align_root_no_ext-all-scores-$partition -nletters $nletters -datafile $data_file -nbestpart $nbestpart -ngroupmax $ngroupmax -o $output_name -ncharbase $min_block_size >/dev/null");
	}

	# Create a Nexus file alignments for each partition
	write_partitions({'PARTITIONS' => $total, 
	                  'NCHAR'      => \%nchar});

	# Clean up unneeded files
	clean_up({'DIRS' => 1});
}

sub write_mdl_command {
	my $settings = shift;

	my @genes = @{$settings->{'GENES'}};

	foreach my $number (0 .. $#genes){

	my $score_file_name = $genes[$number]->{'SCORE_FILE_NAME'};
	my $align_file_name = $genes[$number]->{'ALIGN_FILE_NAME'};
	(my $file_num = $score_file_name) =~ s/.*-(\d+)-\d+/$1/;

	my $start = $genes[$number]->{'GENE_START'} - (($file_num - 1) * $forced_break);
	my $end = $genes[$number]->{'GENE_END'} - (($file_num - 1) * $forced_break);

	my $command;
	$command .= "begin paup;\\nexecute ../$align_file_name;\\n";
	$command .= "set warnroot=no warntree=no warnTsave=no ";
	$command .= "increase=no maxtrees=50 monitor=no notifybeep=no;\\n";
	$command .= "pset gapmode=newstate;\\n" if (!defined($gap_isnt_char));
	$command .= "include $start-$end / only;\\n";
	#$command .= "exclude missambig;\\n" if ($exclude_gaps);
	$command .= "hsearch collapse=no;\\n";
	#$command .= "Pscores 1 / scorefile=$alignment_name-scores-$file_number-".($block_num + 1)." replace=yes;\\n";
	$command .= "Pscores 1 / scorefile=$score_file_name replace=yes;\\n";
	#if ($savetrees) { $command .= "savetrees from=1 to=1 file=$treeFile format=altnexus append=yes;\\n";}

 	#if ($savetrees) {
 	#	$command .= "gettrees file=$treeFile allblocks=yes;\\n";
 	#	$command .= "savetrees file=$alltreeFile format=altnexus replace=yes;\\n";
 	#}
 	$command .= "quit;\\n";
 	$command .= "end;\\n";

	$genes[$number]->{'PAUP_COMMAND'} = $command;
	}
}

sub get_mdl_block_indices {
	my $settings = shift;

	my $file_number = $settings->{'FILE_NUMBER'};
	my $total_genes = $settings->{'TOTAL_GENES'};
	my $align_length = $settings->{'ALIGN_LENGTH'};

	# Determine the start and end indices of all blocks for the given file
	my @blocks;
	foreach my $startblock (1 .. $total_genes) {
	  	foreach my $endblock ($startblock .. $total_genes) {

			my $start_index = (($file_number - 1) * ($forced_break)) + (1 + ($startblock - 1) * $min_block_size);
			my $end_index = (($file_number - 1) * ($forced_break)) + $align_length;
			$end_index = (($file_number - 1) * ($forced_break)) + ($endblock * $min_block_size) if ($endblock < $total_genes);

			push(@blocks, { $start_index => $end_index });
		}
	}
	return @blocks;
}

sub get_num_blocks {
	my $settings = shift;

	my $file_number = $settings->{'FILE_NUMBER'};
	my $total_genes = $settings->{'TOTAL_GENES'};
	my $align_length = $settings->{'ALIGN_LENGTH'};

	# Determine the number of blocks a given file will be broken down into
	my $count = 0;
	foreach my $startblock (1 .. $total_genes) {
	  	foreach my $endblock ($startblock .. $total_genes) {

			my $start_index = (($file_number - 1) * ($forced_break)) + (1 + ($startblock - 1) * $min_block_size);
			my $end_index = (($file_number - 1) * ($forced_break)) + $align_length;
			$end_index = (($file_number - 1) * ($forced_break)) + ($endblock * $min_block_size) if ($endblock < $total_genes);

			$count++;
		}
	}
	return $count;
}

sub get_block_subset {
	my $settings = shift;

	my $file_number = $settings->{'FILE_NUMBER'};

	my @genes = @{$settings->{'GENES'}};
	my $initial_num_genes = scalar(@genes);

	my %break_info = %{$settings->{'BREAK_INFO'}};

	# All job descriptions have been built
	return @genes if (!exists($break_info{$file_number}));

	my $total_genes = $break_info{$file_number}{'TOTAL_GENES'};
	my $align_length = $break_info{$file_number}{'ALIGN_LENGTH'};

	my $size = $settings->{'SIZE'};
	my $genes_needed = $size - scalar(@genes);

	my $start = $settings->{'START'};

	# Determine which forced partition this gene is located in
	my @blocks;
	my $initial_start = 0;
	foreach my $break (sort { $a <=> $b } keys %break_info) {
		last if $break == $file_number;
		$initial_start += $break_info{$break}->{'TOTAL_BLOCKS'};
	}

	# Determine the start and end indices of genes in this partition
	my $count = 0;
	foreach my $startblock (1 .. $total_genes) {
	  	foreach my $endblock ($startblock .. $total_genes) {

			my $start_index = (($file_number - 1) * ($forced_break)) + (1 + ($startblock - 1) * $min_block_size);
			my $end_index = (($file_number - 1) * ($forced_break)) + $align_length;
			$end_index = (($file_number - 1) * ($forced_break)) + ($endblock * $min_block_size) if ($endblock < $total_genes);

			if ($count >= ($start - $initial_start) && $count < (($start - $initial_start) + $genes_needed)) {
				push(@blocks, { $start_index => $end_index });
			}
			$count++;
		}
	}

	my $current_num_genes = scalar(@genes);

	# Create the job description
	foreach my $block_num (0 .. $#blocks) {
		my $block = $blocks[$block_num];

		my $abs_block_num = ($block_num + $start + 1);

		my %info;
		$info{'GENE_START'} = (keys %{$block})[0];
		$info{'GENE_END'} = (values %{$block})[0];
		$info{'PARTITION'} = $count;
		$info{'SCORE_FILE_NAME'} = $align_root_no_ext."-scores-$file_number-$abs_block_num";
		$info{'BLOCK_NUM'} = $abs_block_num;

		$info{'ALIGN_FILE_NAME'} = $align_root_no_ext."-reduced-$file_number.nex";

		push(@genes, \%info);
	}

	# Rerun this method until we have the requested number of job descriptions
	if (scalar(@genes) < $size) {
		my $break = $file_number + 1;
		my $new_size = $size - scalar(@genes);
		my $new_start = $initial_start - ($current_num_genes - $initial_num_genes) + $break_info{$file_number}->{'TOTAL_BLOCKS'};

		return get_block_subset({'FILE_NUMBER'  => $break, 
                                       'START' => $new_start,
                                       'SIZE' => $size,
									   'GENES' => \@genes,
									   'BREAK_INFO' => \%break_info});
	}

	return @genes;
}

sub write_partitions {
	my $settings = shift;

	my $nchar = $settings->{'NCHAR'};
	my $num_partitions = $settings->{'PARTITIONS'};

	my $count = 1;
	my %reduced_partitions;
	foreach my $forced_break (1 .. $num_partitions) {

		# Open MDL partition file
		my $partition_file_name = $partition_dir."$align_root_no_ext-$forced_break"."of$num_partitions.partitions";
		open(my $partition_file, '<', $partition_file_name) 
			or die "Could not open '$partition_file_name': $!.\n";
		chomp(my @data = <$partition_file>);
		close($partition_file);

		my $index;
		foreach my $line (@data) {
			last if ($line =~ /^MDLscore/);
			$index++;
		}
		$index++;

		my @partitions = split(/\s+/, $data[$index]);
		@partitions = splice(@partitions, 2);

		# Doing this makes calculations easier
		$nchar->{0} = $nchar->{1};

		# Determine MDL partitions' starts and ends
		# These indices are for the PIC alignment
		foreach my $index (0 .. $#partitions) {
			my $partition_start = $partitions[$index];

			my $next_partition_start;
			if ($index + 1 < scalar(@partitions)) {
				$next_partition_start = $partitions[$index + 1];
			}
			else {
				$next_partition_start = $nchar->{$forced_break} + 1;
			}

			my $start = $partition_start + (($forced_break - 1) * $nchar->{$forced_break - 1});
			my $end = $next_partition_start + (($forced_break - 1) * $nchar->{$forced_break - 1}) - 1;

			$reduced_partitions{$count}->{'START'} = $start;
			$reduced_partitions{$count}->{'END'} = $end;

			$count++;
		}
	}

	# Determine MDL partitions' starts and ends
	# These indices are for the FULL alignment

	my %stats;
	my %full_partitions;
	foreach my $partition (sort {$a <=> $b} keys %reduced_partitions) {
		my $reduced_start = $reduced_partitions{$partition}->{'START'};
		my $reduced_end = $reduced_partitions{$partition}->{'END'};

		if ($partition == 1) {
			my $next_partition_start = $locations[$reduced_partitions{$partition + 1}->{'START'} - 1];
			my $this_partition_end = $locations[$reduced_partitions{$partition}->{'END'} - 1];
			my $end_offset = floor(($next_partition_start - $this_partition_end - 1) / 2);

			$full_partitions{$partition}->{'START'} = 1;
			$full_partitions{$partition}->{'END'} = $this_partition_end + $end_offset;
		}
		elsif ($partition == scalar(keys %reduced_partitions)) {
			$full_partitions{$partition}->{'START'} = $full_partitions{$partition - 1}->{'END'} + 1;
			$full_partitions{$partition}->{'END'} = length((values %align)[0]); # + 1?
		}
		else {
			my $next_partition_start = $locations[$reduced_partitions{$partition + 1}->{'START'} - 1];
			my $this_partition_end = $locations[$reduced_partitions{$partition}->{'END'} - 1];
			my $end_offset = floor(($next_partition_start - $this_partition_end - 1) / 2);

			$full_partitions{$partition}->{'START'} = $full_partitions{$partition - 1}->{'END'} + 1;
			$full_partitions{$partition}->{'END'} = $this_partition_end + $end_offset;
		}

		$stats{$partition}->{'REDUCED_START'} = $reduced_start;
		$stats{$partition}->{'REDUCED_END'} = $reduced_end;
		$stats{$partition}->{'FULL_START'} = $full_partitions{$partition}->{'START'};
		$stats{$partition}->{'FULL_END'} = $full_partitions{$partition}->{'END'};
		#print "Partition #$partition goes from $reduced_start-$reduced_end when reduced, and $full_partitions{$partition}->{'START'}-$full_partitions{$partition}->{'END'} when expanded.\n";
	}
	print "Alignment has been broken down into ",(scalar(keys %full_partitions))," total partitions.\n";

	# Output potentially useful information on MDL's partitioning
	open(my $stats_file, '>', $align_root_no_ext."-stats.csv");	
	print {$stats_file} "reduced_start, reduced_end, full_start, full_end,\n";
	foreach my $partition (sort { $a <=> $b } keys %stats) {
		print {$stats_file} "$stats{$partition}->{'REDUCED_START'}, $stats{$partition}->{'REDUCED_END'}, $stats{$partition}->{'FULL_START'}, $stats{$partition}->{'FULL_END'}\n";
	}
	close($stats_file);
	print "Gene statistics output to '$align_root_no_ext-stats.csv'.\n";

	# Write the Nexus files which describe the partitioning
	my $free_cpus = get_free_cpus();
	parallelize({'METHOD' => 'write_partition', 
	             'METHOD_ARGS' => {'ALIGN' => \%align, 
				                   'TOKEN' => 'mb', 
								   'PARTITIONS' => \%full_partitions}, 
				 'MAX_FORKS' => $free_cpus, 
				 'TOTAL' => scalar(keys %full_partitions)});

	chdir($gene_dir);	

	# Compress and archive results
	my @gene_file_names = glob("$align_root_no_ext-*-*.nex");
	@gene_file_names = sort { local $a = $a; local $b = $b; $a =~ s/.*-(\d+).nex$/$1/; $b =~ s/.*-(\d+).nex$/$1/; $a <=> $b } @gene_file_names;

	print "Compressing and archiving resulting partitions... ";
	system("tar czf $align_root_no_ext-genes.tar.gz @gene_file_names --remove-files");
	system("mv $align_root_no_ext-genes.tar.gz ..");
	print "done.\n";
}

sub write_partition {
	my $settings = shift;

	my $number = $settings->{'NUM'} + 1;
	my $align = $settings->{'ALIGN'};
	my $partitions = $settings->{'PARTITIONS'};

	my $start = $partitions->{$number}{'START'};
	my $string_start = $start - 1;

	my $end = $partitions->{$number}{'END'} + 1;
	my $string_end = $end - 1;

	my $nchar = $end - $start;

	my $file_name = $gene_dir.$align_root_no_ext."-$start-$string_end.nex";
	unlink($file_name) if (-e $file_name);

	sysopen(my $out, $file_name, O_WRONLY|O_NDELAY|O_CREAT)
		or die "Could not open '$file_name': $!.\n";

	my $line_length = 90000;
	my @taxa = sort { $a cmp $b } (keys %{$align});

	print {$out} "#NEXUS\nbegin data;\n dimensions ntax=".scalar(@taxa)." nchar=$nchar;\n ";
	print {$out} "format datatype=dna interleave=yes gap=- missing=?;\n matrix\n";

	for (my $i = 0; $i < $nchar; $i += $line_length) {
		for (my $j = 0; $j < scalar(@taxa); $j++) {
			my $taxon = $taxa[$j];
			my $interval = $line_length;
			my $current_char = $string_start + $i;

			# Prevents overshoot of alignment output
			
			if ($string_end < $current_char + $interval) {
				$interval = $string_end - $current_char;
			}
			print {$out} "  ".$taxon."\t".substr($align->{$taxon}, $string_start + $i, $interval),"\n";

		}
	}
	print {$out} "\n ;\nend;";
	close($out);

	#print {$out} "begin paup;\ncstatus;\nend;\n";
	#close($out);

	#chomp(my $output = `paup -n $file_name`);
#	chomp(my $output = `paup -n $file_name`);
#
#	(my $pars_inf_chars = $output) =~ s/.*Number of parsimony-informative characters = (\d+).*/$1/s;
#
#	#print $pars_inf_chars," $partition\n";
#	if ($pars_inf_chars == ($reduced_partitions{$number}->{'END'} - $reduced_partitions{$number}->{'START'} + 1)) {
#		#print "expected size and actual size match. ($pars_inf_chars)\n";
#	}
#	else {
#		print "\nPartition $number of ",scalar(keys %reduced_partitions), " expected size and actual size do not match. ($pars_inf_chars)\n";
#		print "",($reduced_partitions{$number}->{'END'} - $reduced_partitions{$number}->{'START'} + 1)," is the expected size.\n";
#	}

}

sub client {
	my ($opt_name, $server_ip) = @_;	

	chdir("/tmp");
	my $paup = "/tmp/paup";

	my $pgrp = $$;

	# Determine this host's IP
	chomp(my $ip = `dig +short myip.opendns.com \@resolver1.opendns.com`); 
	die "Could not establish an IP address for host.\n" if (!defined($ip));

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

	# If the user interrupts the analysis we don't want to leave random files hanging around
	my @unlink;
	$SIG{INT} = sub { unlink(@unlink) };

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
		elsif ($response =~ /NEW: (.*)/) {
			my $args = $1;

			# Parse job information
			(my $paup_command = $args) =~ s/^'.*?' '(.*?)'$/$1/;
			(my $score_file_name = $args) =~ s/^'(.*?)' '.*?'$/$1/;
			
			@unlink = ($score_file_name);

			# Change where alignment is located in command
			$paup_command =~ s/\.\.\///s if ($server_ip ne $ip);

			open(my $std_out, ">&", *STDOUT);
			close(STDOUT);

			# Echo command to PAUP's STDIN
			open(my $paup_pipe, "|-", $paup, "-n"); 
			foreach my $command (split("\\\\n", $paup_command)) {
				print {$paup_pipe} $command,"\n";
			}
			close($paup_pipe);

			open(STDOUT, ">&", $std_out);
			close($std_out);

			# Parse corresponding parsimony score for this block
			open(my $score_file, '<', $score_file_name) 
				or die "Could not open '$score_file_name': $!.\n";
			chomp(my @data = <$score_file>);
			close($score_file);

			(my $score = pop(@data)) =~ s/1\s+(\d+)/$1/;
			unlink(@unlink);

			# Return result to server and ask for a new job
			print {$sock} "RESULT: $score_file_name:$score || NEW: $ip\n";
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

		unlink($0, $paup, glob("*-reduced-*.nex"));
	}
	
	exit(0);
}

sub parallelize {
	my $settings = shift;

	my $method = $settings->{'METHOD'};
	my %method_args = %{$settings->{'METHOD_ARGS'}};
	my $max_forks = $settings->{'MAX_FORKS'};
	my $total = $settings->{'TOTAL'};

	# Methods which can be parallelized by this method
	my %actions = ( 'run_mdl_block'             => \&run_mdl_block,
					'write_partition'           => \&write_partition,
					'write_nexus_file'          => \&write_nexus_file,
					'write_phylip_file_segment' => \&write_phylip_file_segment,
					'write_nexus_file_reduced'  => \&write_nexus_file_reduced);

	my @childPIDs;
	my $running_forks = 0;
	$SIG{CHLD} = 'IGNORE';
	for (my $i = 0; $i < $total; ) {

		# Determine how many processes are currently running
		foreach my $index (reverse(0 .. $#childPIDs)) {

			my $childPID = $childPIDs[$index];
			if (kill 0, $childPID) {
				$running_forks++;
			}
			else {
				splice(@childPIDs, $index, 1);
			}
		}

		# Spawn a new process if running less than maximum
		if ($running_forks < $max_forks) {
			my $pid = fork();
			if ($pid == 0) {
				$method_args{'NUM'} = $i;
				$method_args{'TOTAL'} = $total;
				$actions{$method}->(\%method_args);
				exit(0);
			}
			else {
				push(@childPIDs, $pid);	
				$i++;
			}
		}
		$running_forks = 0;
	}
	foreach my $childPID (@childPIDs) {
		waitpid($childPID, 0);
	}
}

sub parallelize_with_return {
	my $settings = shift;

	my $method = $settings->{'METHOD'};
	my %method_args = %{$settings->{'METHOD_ARGS'}};
	my $max_forks = $settings->{'MAX_FORKS'};
	my $total = $settings->{'TOTAL'};

	# Methods which can be parallelized by this method
	my %actions = ( 'paup_cstatus' => \&paup_cstatus,
	                'tnt_info'     => \&tnt_info,
					'get_informative' => \&get_informative);

	my $select = new IO::Select; 

	my @childPIDs;
	my $running_forks = 0;
	$SIG{CHLD} = 'IGNORE';
	for (my $i = 0; $i < $total; ) {

		# Determine how many processes are currently running
		foreach my $index (reverse(0 .. $#childPIDs)) {
			my $childPID = $childPIDs[$index];
			if (kill 0, $childPID) {
				$running_forks++;
			}
			else {
				splice(@childPIDs, $index, 1);
			}
		}

		# Spawn a new process if running less than maximum
		if ($running_forks < $max_forks) {
			my $pipe = new IO::Pipe;
			my $pid = fork();

			if ($pid == 0) {
				my $TO_PARENT = $pipe->writer();

				$method_args{'NUM'} = $i;
				$method_args{'TOTAL'} = $total;
				$method_args{'TO_PARENT'} = \$TO_PARENT;
				$actions{$method}->(\%method_args);

				exit(0);
			}
			else {
				my $FROM_CHILD = $pipe->reader();
				$select->add($FROM_CHILD);
				push(@childPIDs, $pid);	
				$i++;
			}
		}
		$running_forks = 0;
	}

	my @return;
	while (my @handles = $select->can_read()) {
		my @handles = $select->can_read();
		my $handle = shift(@handles);

		chomp(my @child_response = <$handle>);
		foreach my $response (@child_response) {
			push(@return, $response);
		}
		$select->remove($handle);
		$handle->close;
	}

	return @return;
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

#sub get_free_cpus {
#
#	if ($no_forks) {
#		return 1; # assume that at least one cpu is free
#	}
#	else {
#
#		# Returns a two-member array containing cpu usage observed by the program top,
#		# command is run twice as top's first output is usually inaccurate
#		chomp(my @percent_free_cpu = `top -bn2d0.05 | grep "Cpu(s)"`);
#
#		my $percent_free_cpu = pop(@percent_free_cpu);
#		my $test = $percent_free_cpu;
#		#$percent_free_cpu =~ s/.*?(\d+\.\d)%ni,\s*(\d+\.\d)%id.*/$1 + $2/; # also includes %nice as free 
#		$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
#		$percent_free_cpu = eval($percent_free_cpu);
#
#		my $total_cpus = `grep 'cpu' /proc/stat | wc -l` - 1;
#		die "$test\n" if (!defined($percent_free_cpu));
#
#		my $free_cpus = ceil($total_cpus * $percent_free_cpu / 100);
#
#		if ($free_cpus == 0) {
#			$free_cpus = 1; # assume that at least one cpu can be used
#		}
#		
#		return $free_cpus;
#	}
#}

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

sub INT_handler {
	clean_up({'DIRS' => 1});
	exit(0);
}

sub clean_up {
	my $settings = shift;

	my $remove_dirs = $settings->{'DIRS'};
	my $current_dir = getcwd();

	chdir($project_name);
#	chdir($alignment_root);
#	unlink(glob($gene_dir."$alignment_name*"));
	unlink(glob($score_dir."$align_root_no_ext*"));
#	unlink(glob($phylip_dir."$alignment_name*"));
#	#unlink(glob($partition_dir."$alignment_name*"));
#	unlink(glob($alignment_name."-unreduced-*"));

	if ($remove_dirs) {
		rmdir($gene_dir);
		rmdir($score_dir);
		rmdir($phylip_dir);
		#rmdir($partition_dir);
	}
	chdir($current_dir);
}

sub usage {
	return "Usage: alignment-breakdown.pl [ALIGNMENT] [-b MINIMUM BLOCK LENGTH]\n";
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

sub logger {
	my $msg = shift;

	my $time = "[".localtime(time())."]";

	# Allow for new lines and carriage returns before message
	if ($msg =~ s/^\n//) {
		$time = "\n$time";
	}
	elsif ($msg =~ s/^\r//) {
		$time = "\r$time";
	}

	print "$time $msg\n"; 
}

sub help {
	return "alignment-breakdown.pl help\n";
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
	die if ($client);

	die "Could not find the following executable: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path));
	return $exec_path;
}
