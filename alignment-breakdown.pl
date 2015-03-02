#!/usr/bin/perl
use strict;
use warnings;
#use Cwd;
#use Cwd 'abs_path';
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
my $command_dir   = "mdl-commands/";
my $partition_dir = "mdl-partitions/";

mkdir($gene_dir)      or die "Could not create '$gene_dir': $!.\n"      if (!-e $gene_dir);
mkdir($score_dir)     or die "Could not create '$score_dir': $!.\n"     if (!-e $score_dir);
mkdir($phylip_dir)    or die "Could not create '$phylip_dir': $!.\n"    if (!-e $phylip_dir);
mkdir($command_dir)   or die "Could not create '$command_dir': $!.\n"   if (!-e $command_dir);
mkdir($partition_dir) or die "Could not create '$partition_dir': $!.\n" if (!-e $partition_dir);

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

	# Run PAUP via a pipe to save time
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
##	parallelize({'METHOD'      => 'write_nexus_file_reduced', 
##	             'METHOD_ARGS' => {'ALIGN' => \%align, 
##				                   'TOKEN' => 'reduced'}, 
##				 'MAX_FORKS'   => $free_cpus, 
##				 'TOTAL'       => $total_output_files});

	print "\n", $num_pars_inf_chars, " total parsimony-informative sites found for '$align'.\n\n";

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
	chomp(my $server_ip = `wget -qO- ipecho.net/plain`); # Returns the external IP address of this computer
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

			# Send reduced alignments to the machine
			if ($machine ne $server_hostname) {
				system("scp -q *-reduced-* $machine:/tmp");
			}

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

								# Determine which reduced output file we are in
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
								
								my $time = time();
								my $free_cpus = get_free_cpus();
								parallelize({'METHOD'      => 'write_mdl_command', 
											 'METHOD_ARGS' => {'GENES' => \@gene_subset}, 
											 'MAX_FORKS'   => $free_cpus, 
											 'TOTAL'       => scalar(@gene_subset)});

								print "done. ($start - ".($start + scalar(@gene_subset) - 1).")\n";
								print "  Total execution time: ", secs_to_readable({'TIME' => time() - $time}), "\n";

								$subset_count++;
							}
							
							my $job = shift(@gene_subset);
							#my $phylip_file_name = $job->{'PHY_FILE_NAME'};
							my $score_file_name = $job->{'SCORE_FILE_NAME'};
							my $com_file_name = $job->{'COM_FILE_NAME'};
							my $align_file_name = $job->{'ALIGN_FILE_NAME'};

							# TODO: send actual commands instead of file to save time
							# Check whether the client is remote or local, send it needed files if remote
							if ($client_ip ne $server_ip) {

								$SIG{CHLD} = 'IGNORE';

								# Fork to perform the file transfer and prevent stalling the server
								my $pid = fork(); 
								if ($pid == 0) {
									
									chdir($command_dir);
									send_file({'FILE_PATH' => $com_file_name, 'FILE_HANDLE' => $client});						

									print {$client} "NEW: -n '$com_file_name'\n";
									exit(0);
								}
							}
							else {
								print {$client} "CHDIR: ".abs_path($command_dir)."\n";
								print {$client} "NEW: -n '$com_file_name'\n";
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

#		my $time = time();
#		my $free_cpus = get_free_cpus();
#		my @genes = get_submit_file_contents_flat({'GENES' => \%genes});
#		$time = time();
#		my $free_cpus = get_free_cpus();
#		parallelize({'METHOD'      => 'run_mdl_block', 
#					 'METHOD_ARGS' => {'ALIGN' => \%align, 
#									   'GENES' => \@genes}, 
#									   #'GENES' => \%genes}, 
#					 'MAX_FORKS'   => $free_cpus, 
#					 #'TOTAL'       => 1});
#					 'TOTAL'       => scalar(@genes)});
#		print "All submitted jobs are complete.\n\n";
#		print "Total execution time: ", secs_to_readable({'TIME' => time() - $time}), "\n\n";
		#die;
	#}

	#chdir($alignment_root);
	chdir($project_name);
	#parse_input({'GET_ALIGN' => 1});
	#get_informative_chars({'ALIGN' => \%align});

	#unlink(glob("$alignment_name-reduced*"));

	my $ntax = scalar(keys %align);

	#my $mdl = check_path_for_exec("mdl");
	
	# Run each partition defined from the forced break points (if any) through MDL
	foreach my $partition (1 .. $total) {
		#my $data_file = "$alignment_name-reduced-$partition"."of$total.nex";
		my $data_file = "$align_root_no_ext-reduced-$partition"."of$total.nex";
		#my $output_name = $partition_dir."$alignment_name-$partition"."of$total.partitions";
		my $output_name = $partition_dir."$align_root_no_ext-$partition"."of$total.partitions";
		#system("./mdl -ntax $ntax -nchar $nchar{$partition} -scorefile mdl-scores/$alignment_name-all-scores-$partition -nletters $nletters -datafile $data_file -nbestpart $nbestpart -ngroupmax $ngroupmax -o $output_name -ncharbase $min_block_size >/dev/null");
		#system("$mdl -ntax $ntax -nchar $nchar{$partition} -scorefile mdl-scores/$alignment_name-all-scores-$partition -nletters $nletters -datafile $data_file -nbestpart $nbestpart -ngroupmax $ngroupmax -o $output_name -ncharbase $min_block_size >/dev/null");
		system("$mdl -ntax $ntax -nchar $nchar{$partition} -scorefile mdl-scores/$align_root_no_ext-all-scores-$partition -nletters $nletters -datafile $data_file -nbestpart $nbestpart -ngroupmax $ngroupmax -o $output_name -ncharbase $min_block_size >/dev/null");
	}

	write_partitions({'PARTITIONS' => $total, 
	                  'NCHAR'      => \%nchar});

	clean_up({'DIRS' => 1});
}

#sub run_mdl_block {
#	my $settings = shift;
#
#	my $number = $settings->{'NUM'};
#	#my %genes = %{$settings->{'GENES'}};
#	my @genes = @{$settings->{'GENES'}};
#	
#	#$number = 27630;
#
#	my $phylip_file_name = "../mdl-phylip/".$genes[$number]->{'PHY_FILE_NAME'};
#	my $score_file_name = $genes[$number]->{'SCORE_FILE_NAME'};
#
#
#	write_phylip_file_segment({'NUM'   => $number,
#							   'ALIGN' => \%align, 
#						       'GENES' => \@genes}); 
#						       #'GENES' => \%genes}); 
#
#	chdir($score_dir);
#
#	my $seed = int(rand(100000));
#	system("../parsimonator-SSE3 -s '$phylip_file_name' -n '$score_file_name' -p $seed >/dev/null");
#
#	
##	unlink($phylip_file_name);
##	unlink("RAxML_parsimonyTree.".$score_file_name.".0");
#	#unlink($phylip_file_name) if $number != 27630;
#	#unlink("RAxML_parsimonyTree.".$score_file_name.".0") if $number != 27630;
#
#	(my $partition = $score_file_name) =~ s/.*-scores-(\d+)-\d+.*/$1/;
#
##	open(my $joined_score_file, '>>', $alignment_name."-all-scores-$partition") 
##		or die "Could not open '".$alignment_name."-all-scores-$partition': $!.\n";
#	open(my $joined_score_file, '>>', $align_root_no_ext."-all-scores-$partition") 
#		or die "Could not open '".$align_root_no_ext."-all-scores-$partition': $!.\n";
#	flock($joined_score_file, LOCK_EX) or die "Could not lock '$align_root_no_ext-all-scores-$partition': $!.\n";
#	seek($joined_score_file, 0, SEEK_END);
#
#	open(my $score_file, '<', "RAxML_info.".$score_file_name) 
#		or die "Could not open 'RAxML_info.$score_file_name': $!.\n";
#	chomp(my @data = <$score_file>);
#	close($score_file);
#
#	(my $score = pop(@data)) =~ s/.*with length (\d+) computed.*/$1/;
#	(my $tree_number = $score_file_name) =~ s/.*-scores-\d+-(\d+).*/$1/;
#
#	print {$joined_score_file} "$tree_number : $score\n";
#
#	close($joined_score_file);
#
#	#unlink("RAxML_info.".$score_file_name);
#	#unlink("RAxML_info.".$score_file_name) if $number != 27630;
#
#	$number++;
#	if ($number % 1000 == 0 || $number == scalar(@genes)) {
#		my $num_digits = get_num_digits({'NUMBER' => scalar(@genes)});
#		printf("  Analyses complete: %".$num_digits."d/%d.\n", $number, scalar(@genes));
#	}
#
#}

sub write_mdl_command {
	my $settings = shift;

	my $number = $settings->{'NUM'};
	my @genes = @{$settings->{'GENES'}};

	#(my $gene_number = $genes[$number]->{'SCORE_FILE_NAME'}) =~ s/\Q$alignment_name\E-scores-//;
	my $score_file_name = $genes[$number]->{'SCORE_FILE_NAME'};
	my $align_file_name = $genes[$number]->{'ALIGN_FILE_NAME'};
	(my $file_num = $score_file_name) =~ s/.*-(\d+)-\d+/$1/;
	
#	my $start = $genes[$number]->{'GENE_START'};
#	my $end = $genes[$number]->{'GENE_END'};
	my $start = $genes[$number]->{'GENE_START'} - (($file_num - 1) * $forced_break);
	my $end = $genes[$number]->{'GENE_END'} - (($file_num - 1) * $forced_break);

	#my $file_name = $phylip_dir.$genes[$number]->{'PHY_FILE_NAME'};
	my $file_name = $command_dir.$genes[$number]->{'COM_FILE_NAME'};
	unlink($file_name) if (-e $file_name);

	sysopen(my $out, $file_name, O_WRONLY|O_NDELAY|O_CREAT)
		or die "Could not open '$file_name': $!.\n";

	my $command;
	#$command .= "begin paup;\nexecute $align_file_name\n";
	$command .= "begin paup;\nexecute ../$align_file_name;\n";
	$command .= "set warnroot=no warntree=no warnTsave=no ";
	$command .= "increase=no maxtrees=50 monitor=no notifybeep=no;\n";
	$command .= "pset gapmode=newstate;\n" if (!defined($gap_isnt_char));
	$command .= "include $start-$end / only;\n";
	#$command .= "exclude missambig;\n" if ($exclude_gaps);
	$command .= "hsearch collapse=no;\n";
	#$command .= "Pscores 1 / scorefile=$alignment_name-scores-$file_number-".($block_num + 1)." replace=yes;\n";
	#$command .= "Pscores 1 / scorefile=$score_file_name replace=yes;\n";
	#$command .= "Pscores 1 / scorefile=../mdl-scores/$score_file_name replace=yes;\n";
	$command .= "Pscores 1 / scorefile=$score_file_name replace=yes;\n";
	#if ($savetrees) { $command .= "savetrees from=1 to=1 file=$treeFile format=altnexus append=yes;\\n";}

 	#if ($savetrees) {
 	#	$command .= "gettrees file=$treeFile allblocks=yes;\\n";
 	#	$command .= "savetrees file=$alltreeFile format=altnexus replace=yes;\\n";
 	#}
 	$command .= "quit;\n";
 	$command .= "end;\n";

	print {$out} $command;
	close($out);
}

sub get_mdl_block_indices {
	my $settings = shift;

	my $file_number = $settings->{'FILE_NUMBER'};
	my $total_genes = $settings->{'TOTAL_GENES'};
	my $align_length = $settings->{'ALIGN_LENGTH'};

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

	return @genes if (!exists($break_info{$file_number}));
	#return if (!exists($break_info{$file_number}));

	my $total_genes = $break_info{$file_number}{'TOTAL_GENES'};
	my $align_length = $break_info{$file_number}{'ALIGN_LENGTH'};

	my $size = $settings->{'SIZE'};
	my $genes_needed = $size - scalar(@genes);

	my $start = $settings->{'START'};

	#print "\ncurrent genes: ".scalar(@genes).", start: $start, size: $size, genes_needed: $genes_needed, block: $file_number\n";

	my @blocks;
	#my $count = 0;
	my $initial_start = 0;
	foreach my $break (sort { $a <=> $b } keys %break_info) {
		last if $break == $file_number;
		$initial_start += $break_info{$break}->{'TOTAL_BLOCKS'};
	}
	#print $count,"\n";

	my $count = 0;
	foreach my $startblock (1 .. $total_genes) {
	  	foreach my $endblock ($startblock .. $total_genes) {

			my $start_index = (($file_number - 1) * ($forced_break)) + (1 + ($startblock - 1) * $min_block_size);
			my $end_index = (($file_number - 1) * ($forced_break)) + $align_length;
			$end_index = (($file_number - 1) * ($forced_break)) + ($endblock * $min_block_size) if ($endblock < $total_genes);

			#if ($count >= $start && $count < ($start + $size)) {
			#if ($count >= $start && $count < ($start + $genes_needed)) {
			#if ($count >= ($start - $initial_start) && $count < ($start + $genes_needed)) {
			if ($count >= ($start - $initial_start) && $count < (($start - $initial_start) + $genes_needed)) {
				push(@blocks, { $start_index => $end_index });
			}
			$count++;
		}
	}

	my $current_num_genes = scalar(@genes);
	#my $current_gene_total = $#genes;

	#my $initial_start = $count - scalar(@blocks);

	foreach my $block_num (0 .. $#blocks) {
		my $block = $blocks[$block_num];

		my $abs_block_num = ($block_num + $start + 1);
		#my $phylip_file_name = $alignment_name."-$file_number-$abs_block_num.phy";
		my $phylip_file_name = $align_root_no_ext."-$file_number-$abs_block_num.phy";

		my %info;
		$info{'PHY_FILE_NAME'} = $phylip_file_name;
		$info{'GENE_START'} = (keys %{$block})[0];
		$info{'GENE_END'} = (values %{$block})[0];
		$info{'PARTITION'} = $count;
		($info{'SCORE_FILE_NAME'} = $info{'PHY_FILE_NAME'}) =~ s/-\Q$file_number\E/-scores-$file_number/;
		$info{'SCORE_FILE_NAME'} =~ s/\.phy$//;
		#$info{'SCORE_FILE_NAME'} =~ s/\.phy$//;
		$info{'BLOCK_NUM'} = $abs_block_num;

		#$info{'ALIGN_FILE_NAME'} = $alignment_name."-reduced-$file_number.nex";
		$info{'ALIGN_FILE_NAME'} = $align_root_no_ext."-reduced-$file_number.nex";
		($info{'COM_FILE_NAME'} = $info{'SCORE_FILE_NAME'}) =~ s/-scores-/-commands-/;

		push(@genes, \%info);
	}


	if (scalar(@genes) < $size) {
		#print "more genes needed.\n";
		my $break = $file_number + 1;
		my $new_size = $size - scalar(@genes);
		my $new_start = $initial_start - ($current_num_genes - $initial_num_genes) + $break_info{$file_number}->{'TOTAL_BLOCKS'};
		return get_block_subset({'FILE_NUMBER'  => $break, 
		#get_block_subset({'FILE_NUMBER'  => $break, 
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

		#my $partition_file_name = $partition_dir."$alignment_name-$forced_break"."of$num_partitions.partitions";
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

		$nchar->{0} = $nchar->{1};

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
		#$stats{$partition}->{'REDUCED_LENGTH'} = $reduced_end - $reduced_start;
		$stats{$partition}->{'FULL_START'} = $full_partitions{$partition}->{'START'};
		$stats{$partition}->{'FULL_END'} = $full_partitions{$partition}->{'END'};
		#$stats{$partition}->{'FULL_LENGTH'} = $full_partitions{$partition}->{'END'} - $full_partitions{$partition}->{'START'};
		#$stats{$partition}->{'PARS_FREQ'} = $stats{$partition}->{'REDUCED_LENGTH'} / $stats{$partition}->{'FULL_LENGTH'};

		#print "Partition #$partition goes from $reduced_start-$reduced_end when reduced, and $full_partitions{$partition}->{'START'}-$full_partitions{$partition}->{'END'} when expanded.\n";
	}
	print "Alignment has been broken down into ",(scalar(keys %full_partitions))," total partitions.\n";

	#open(my $stats_file, '>', $alignment_name."-stats.csv");	
	open(my $stats_file, '>', $align_root_no_ext."-stats.csv");	
	#print "reduced_start, reduced_end, reduced_length, full_start, full_end, full_length, pars_freq\n";
	print {$stats_file} "reduced_start, reduced_end, full_start, full_end,\n";
	foreach my $partition (sort { $a <=> $b } keys %stats) {
		#print {$stats_file} "$stats{$partition}->{'REDUCED_START'}, $stats{$partition}->{'REDUCED_END'}, $stats{$partition}->{'REDUCED_LENGTH'}, $stats{$partition}->{'FULL_START'}, $stats{$partition}->{'FULL_END'}, $stats{$partition}->{'FULL_LENGTH'}, $stats{$partition}->{'PARS_FREQ'}\n";
		print {$stats_file} "$stats{$partition}->{'REDUCED_START'}, $stats{$partition}->{'REDUCED_END'}, $stats{$partition}->{'FULL_START'}, $stats{$partition}->{'FULL_END'}\n";
	}
	close($stats_file);
	#print "Gene statistics output to '$alignment_name-stats.csv'.\n";
	print "Gene statistics output to '$align_root_no_ext-stats.csv'.\n";

	my $free_cpus = get_free_cpus();
	parallelize({'METHOD' => 'write_partition', 
	             'METHOD_ARGS' => {'ALIGN' => \%align, 
				                   'TOKEN' => 'mb', 
								   'PARTITIONS' => \%full_partitions}, 
				 'MAX_FORKS' => $free_cpus, 
				 'TOTAL' => scalar(keys %full_partitions)});

	chdir($gene_dir);	

	#my @gene_file_names = glob("$alignment_name-*-*.nex");
	my @gene_file_names = glob("$align_root_no_ext-*-*.nex");
	@gene_file_names = sort { local $a = $a; local $b = $b; $a =~ s/.*-(\d+).nex$/$1/; $b =~ s/.*-(\d+).nex$/$1/; $a <=> $b } @gene_file_names;

	print "Compressing and archiving resulting partitions... ";
	#system("tar czf $alignment_name-genes.tar.gz @gene_file_names --remove-files");
	system("tar czf $align_root_no_ext-genes.tar.gz @gene_file_names --remove-files");
	#system("mv $alignment_name-genes.tar.gz $alignment_root");
	#system("mv $alignment_name-genes.tar.gz ..");
	system("mv $align_root_no_ext-genes.tar.gz ..");
	print "done.\n";

}

sub write_partition {
	my $settings = shift;

	my $number = $settings->{'NUM'} + 1;
	my $align = $settings->{'ALIGN'};
	#my $total = $settings->{'TOTAL'};
	#my $token = $settings->{'TOKEN'};
	my $partitions = $settings->{'PARTITIONS'};

	my $start = $partitions->{$number}{'START'};
	my $string_start = $start - 1;

	my $end = $partitions->{$number}{'END'} + 1;
	my $string_end = $end - 1;

	my $nchar = $end - $start;

	#my $file_name = $gene_dir.$alignment_name."-$start-$string_end.nex";
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
	#my $server_ip = shift;	
	my ($opt_name, $server_ip) = @_;	

	#chdir($ENV{'HOME'});
	chdir("/tmp");
	
	my $paup = "/tmp/paup";

	my $pgrp = $$;
	#setpgrp($$, $pgrp);
	#setpgrp(0, $pgrp);
	#$SIG{CHLD} = 'IGNORE';

	chomp(my $ip = `wget -qO- ipecho.net/plain`); 
	die "Could not establish an IP address for host.\n" if (not defined $ip);

	# Spawn more clients

	my @pids;
	my $total_forks = get_free_cpus();
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

	# If the user interrupts the analysis we don't want to leave random files in their home directory
	my @unlink;
	$SIG{INT} = sub { unlink(@unlink) };

#	my $sock = new IO::Socket::INET(
#		PeerAddr  => $server_ip.":".$port,
#		Proto => 'tcp') 
#	or die "Could not connect to server socket: $!.\n"; 
	my $sock = new IO::Socket::INET(
		PeerAddr  => $server_ip.":".$port,
		Proto => 'tcp') 
	or exit(0); 
	$sock->autoflush(1);


	#sleep(int(rand(10)));

	#print "Client succesfully generated with IP: $ip.\n";

	print {$sock} "NEW: $ip\n";
	while (chomp(my $response = <$sock>)) {
		#chomp($response);

		#if ($response =~ /FILE: (.*)/) {
		if ($response =~ /SEND_FILE: (.*)/) {
			my $file_name = $1;
			#receive_file($file_name, $sock);
			receive_file({'FILE_PATH' => $file_name, 'FILE_HANDLE' => $sock});						
			#receive_file({'FILE_PATH' => $file_name, 'FILE_HANDLE' => \$sock});						
		}
		elsif ($response =~ /CHDIR: (.*)/) {
			chdir($1);
		}
		elsif ($response =~ /NEW: (.*)/) {
			my $args = $1;

		#	chomp(my $host = `hostname`);
		#	print "$host : $fork_id\n";

#			(my $score_file_name = $args) =~ s/.*-n '(.*?)'.*/$1/;
#			(my $phylip_file_name = $args) =~ s/.*-s '(.*?)'.*/$1/;
#			(my $seed = $args) =~ s/.*-p (\d+).*/$1/;
			
			(my $com_file_name = $args) =~ s/.*-n '(.*?)'.*/$1/;
			(my $score_file_name = $com_file_name) =~ s/-commands-/-scores-/;

			#@unlink = ($phylip_file_name, "RAxML_info.".$score_file_name, "RAxML_parsimonyTree.".$score_file_name.".0");
			@unlink = ($com_file_name, $score_file_name);

			# Change where alignment is located in command file
			if ($server_ip ne $ip) {

				open(my $com_file, "<", $com_file_name);
				my @data = <$com_file>;
				close($com_file);

				foreach my $line (@data) {
					$line =~ s/\.\.\///;
				}

				open($com_file, ">", $com_file_name);
				print {$com_file} @data;
				close($com_file);
			}

			#if (cwd() eq $ENV{HOME}) {
#			if ($cwd =~ /\Q$ENV{HOME}\E$/) {
#				system("./parsimonator-SSE3 $args >/dev/null");
#			}
#			else {
##				system("../parsimonator-SSE3 $args >/dev/null") == 0 or die "\nCall to parsimonator failed: $!. ($args).\n";
#	#			open(my $std_out, ">&", *STDOUT);
#	#			close(STDOUT);
#				
#				#(my $args = $args) =~ s/'//g;
#
#				#open(my $parsimonator_pipe, "-|", "../parsimonator-SSE3", $args); 
#				#system("../parsimonator-SSE3", split(/\s+/, $args));
#				#system("../parsimonator-SSE3", "-n", $score_file_name, "-s", $phylip_file_name, "-p", $seed);
#				#system("../paup", "-n", $com_file_name);
#				#close($parsimonator_pipe);
#
#	#			open(STDOUT, ">&", $std_out);
#	#			close($std_out);
#			}

			open(my $std_out, ">&", *STDOUT);
			close(STDOUT);

			system($paup, "-n", $com_file_name);

			open(STDOUT, ">&", $std_out);
			close($std_out);

			(my $partition = $com_file_name) =~ s/.*-commands-(\d+)-\d+.*/$1/;

			open(my $score_file, '<', $score_file_name) 
				or die "Could not open '$score_file_name': $!.\n";
			chomp(my @data = <$score_file>);
			close($score_file);

			(my $score = pop(@data)) =~ s/1\s+(\d+)/$1/;
			unlink(@unlink);

			print {$sock} "RESULT: $score_file_name:$score || NEW: $ip\n";
		}
		elsif ($response eq "HANGUP") {
			last;
		}
	}

	if ($$ == $pgrp) {
		foreach my $pid (@pids) {
			waitpid($pid, 0);
		}

		unlink($0, $paup, glob("*-reduced-*.nex"));

		#print "Exiting: '",last_in_pgrp({'PGRP' => $pgrp, 'PID' => $$}),"'\n";
		#print `ps --ppid $pgrp`;
		#cleanup
	}
	
	exit(0);
}

sub parallelize {
	my $settings = shift;

	my $method = $settings->{'METHOD'};
	my %method_args = %{$settings->{'METHOD_ARGS'}};
	my $max_forks = $settings->{'MAX_FORKS'};
	my $total = $settings->{'TOTAL'};

	#print "Running '$method' $total total times with $max_forks forks.\n";

	my %actions = ( 'run_mdl_block'             => \&run_mdl_block,
					'write_partition'           => \&write_partition,
					'write_nexus_file'          => \&write_nexus_file,
					'write_mdl_command'         => \&write_mdl_command,
					'write_phylip_file_segment' => \&write_phylip_file_segment,
					'write_nexus_file_reduced'  => \&write_nexus_file_reduced);

	#use Time::HiRes qw(time);

	my $time = time();

	my @childPIDs;
	my $running_forks;
	$SIG{CHLD} = 'IGNORE';
	for (my $i = 0; $i < $total; ) {
		foreach my $childPID (@childPIDs) {
			$running_forks++ if (kill 0, $childPID);
		}
		if (!defined $running_forks || $running_forks < $max_forks) {
			my $time = time();
			my $pid = fork();
	#		print "secs to fork: ",(time() - $time),"\n";
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
	#	if ($method eq 'write_phylip_file_segment' && $i == 10) {
	#		print "total secs: ",(time() - $time),"\n";
	#		die;
	#	}
		undef($running_forks);
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

	my %actions = ( 'paup_cstatus' => \&paup_cstatus,
	                'tnt_info'     => \&tnt_info,
					'get_informative' => \&get_informative);

	my $select = new IO::Select; 

	my @childPIDs;
	my $running_forks;
	$SIG{CHLD} = 'IGNORE';
	for (my $i = 0; $i < $total; ) {

		foreach my $childPID (@childPIDs) {
			$running_forks++ if (kill 0, $childPID);
		}
		if (!defined $running_forks || $running_forks < $max_forks) {
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
		undef($running_forks);
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

	if ($no_forks) {
		return 1; # assume that at least one cpu is free
	}
	else {

		# Returns a two-member array containing cpu usage observed by the program top,
		# command is run twice as top's first output is usually inaccurate
		chomp(my @percent_free_cpu = `top -bn2d0.05 | grep "Cpu(s)"`);

		my $percent_free_cpu = pop(@percent_free_cpu);
		my $test = $percent_free_cpu;
		#$percent_free_cpu =~ s/.*?(\d+\.\d)%ni,\s*(\d+\.\d)%id.*/$1 + $2/; # also includes %nice as free 
		$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
		$percent_free_cpu = eval($percent_free_cpu);

		my $total_cpus = `grep 'cpu' /proc/stat | wc -l` - 1;
		die "$test\n" if (!defined($percent_free_cpu));

		my $free_cpus = ceil($total_cpus * $percent_free_cpu / 100);

		if ($free_cpus == 0) {
			$free_cpus = 1; # assume that at least one cpu can be used
		}
		
		return $free_cpus;
	}
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
	unlink(glob($command_dir."$align_root_no_ext*"));
#	#unlink(glob($partition_dir."$alignment_name*"));
#	unlink(glob($alignment_name."-unreduced-*"));

	if ($remove_dirs) {
		rmdir($gene_dir);
		rmdir($score_dir);
		rmdir($phylip_dir);
		rmdir($command_dir);
		#rmdir($partition_dir);
	}
	chdir($current_dir);
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

	#my $file_path = $settings->{'FILE_PATH'};
	#my $file_handle = $settings->{'FILE_HANDLE'};
	my $file_path = $settings->{'FILE_PATH'};
	my $file_handle = $settings->{'FILE_HANDLE'};

	my $hash = hashsum({'FILE_PATH' => $file_path});
	print {$file_handle} "SEND_FILE: $file_path\n";
	#print {$$file_handle} "SEND_FILE: $file_path\n";

	open(my $file, "<", $file_path) or die "Couldn't open file '$file_path': $!.\n";
	while (<$file>) {
		print {$file_handle} $_;
		#print {$$file_handle} $_;
	}
	#print {$handle} $_ while (local $_ = <$file>);
	close($file);

	#print {$file_handle} "END_FILE: $hash\n";
	print {$file_handle} " END_FILE: $hash\n";
	#print {$$file_handle} "END_FILE: $hash\n";
	#print {$$file_handle} " END_FILE: $hash\n";
}

sub receive_file {
	my $settings = shift;

	#my $file_path = $settings->{'FILE_PATH'};
	#my $file_handle = $settings->{'FILE_HANDLE'};
	my $file_path = $settings->{'FILE_PATH'};
	my $file_handle = $settings->{'FILE_HANDLE'};

	my $check_hash;
	open(my $file, ">", $file_path);
	while (<$file_handle>) {
	#while (<$$file_handle>) {
		#if ($_ =~ /(.*)END_FILE: (\S+)/) {
		if ($_ =~ /(.*) END_FILE: (\S+)/) {
			#print {$file_handle} $1;
			#print {$$file_handle} $1;
			print {$file} $1;
			$check_hash = $2;
			last;
		}
		else {
			#print {$file_handle} $_;
			#print {$$file_handle} $_;
			print {$file} $_;
		}
	}
	close($file);
	#close($$file_handle);

	#my $hash = hashsum($file_path);
	my $hash = hashsum({'FILE_PATH' => $file_path});
	if ($hash ne $check_hash) {
		print "Unsuccessful file transfer, checksums do not match.\n'$hash' - '$check_hash'\n"; # hopefully this never pops up
	}
}

sub usage {
	return "Usage: alignment-breakdown.pl [ALIGNMENT] [-b MINIMUM BLOCK LENGTH\n";
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
