#!/usr/bin/perl
use strict;
use warnings;
use Fcntl;
use Fcntl qw(:flock);
use POSIX;
use IO::Pipe;
use IO::Handle;
use IO::Select;
use IO::Socket;
use Getopt::Long;
use Cwd qw(abs_path);
use Time::HiRes qw(time usleep);

my $os_name = $^O;

# prevents printing a recursion warning which appear with large datasets
no warnings 'recursion';

# Turn on autoflush
$|++;

# Maximum number of CPUs to use
my $max_forks;

# Server port
my $port = 10001;

# Stores executing machine hostnames
my @machines;
my %machines;

# Path to text file containing computers to run on
my $machine_file_path;

# Where this script is located 
my $script_path = abs_path($0);

# MDL settings
my $gap_as_char;
#my $gap_isnt_char;
my $exclude_gaps;
my $nletters = 4;
my $nbestpart = 1;
my $ngroupmax = 10000;

# General script settings
my $forced_break;
my $min_block_size;

# Stores total number of parsimony-informative characters found by script
my $num_pars_inf_chars;

# How the script was called
my $invocation = "perl mdl.pl @ARGV";

# Name of output directory
my $project_name = "mdl-".int(time());
#my $project_name = "mdl-dir";

# Read commandline settings
GetOptions(
	"block-size|b=i"   => \$min_block_size,
	"nletters|l=i"     => \$nletters,
	#"nbestpart|p=i"    => \$nbestpart,
	#"ngroupmax|m=i"    => \$ngroupmax,
	"gap-as-char|g"    => \$gap_as_char,
	#"gap-isnt-char|g"  => \$gap_isnt_char,
    "n-threads|T=i"    => \$max_forks,
    "out-dir|o=s"      => \$project_name,
	"exclude-gaps|e"   => \$exclude_gaps,
	"forced-break|f=i" => \$forced_break,
	"machine-file=s"   => \$machine_file_path,
	"port|p=i"         => \$port,
	"server-ip=s"      => \&client, # for internal usage only
	"help|h"           => sub { print &help; exit(0); },
	"usage"            => sub { print &usage; exit(0); },
);

# Get paths to required executables
my $mdl = check_path_for_exec("mdl");
my $paup = check_path_for_exec("paup");

my $align = shift(@ARGV);

# Some error checking
die "You must specify an alignment file.\n\n", &usage if (!defined($align));
die "Could not locate '$align', perhaps you made a typo.\n" if (!-e $align);
die "You must specify a minimum block length (-b).\n\n", &usage if (!defined($min_block_size));
die "The minimum block size must be greater than 0!\n" if ($min_block_size <= 0);
die "Forced breakpoints can't be negative length!\n" if (defined($forced_break) && $forced_break <= 0);
die "Could not locate '$machine_file_path'.\n" if (defined($machine_file_path) && !-e $machine_file_path);
die "The minimum block size can not be greater than or equal to the forced breakpoint length.\n" if (defined($forced_break) && $min_block_size >= $forced_break);

print "WARNING: It is not recommended to have forced breakpoints placed so frequently unless you are analyzing a very large dataset.\n" if (defined($forced_break) && $forced_break <= 100);

# Adjust MDL settings depending on user input
#$gap_isnt_char = 1 if (defined($exclude_gaps));
$gap_as_char = 0 if (defined($exclude_gaps));
#$nletters++      if (!defined($exclude_gaps));
$nletters++      if (!defined($exclude_gaps) && $gap_as_char);

# Extract name information from input file
(my $align_root = $align) =~ s/.*\/(.*)/$1/;
(my $align_root_no_ext = $align) =~ s/(.*\/)?(.*)\..*/$2/;

# Initialize working directory
# Remove conditional eventually
mkdir($project_name) || die "Could not create '$project_name'$!.\n" if (!-e $project_name);

my $align_abs_path = abs_path($align);
# Remove conditional eventually
run_cmd("ln -s $align_abs_path $project_name/$align_root") if (! -e "$project_name/$align_root");

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

chdir($project_name);

# Define and initialize directories
my $gene_dir      = "mdl-genes/";
my $score_dir     = "mdl-scores/";
my $partition_dir = "mdl-partitions/";

mkdir($gene_dir)      or die "Could not create '$gene_dir': $!.\n"      if (!-e $gene_dir);
mkdir($score_dir)     or die "Could not create '$score_dir': $!.\n"     if (!-e $score_dir);
mkdir($partition_dir) or die "Could not create '$partition_dir': $!.\n" if (!-e $partition_dir);

# Change how Ctrl+C is interpreted to allow for clean up
$SIG{'INT'} = 'INT_handler';

# Print the current script settings
print "\nScript was called as follows:\n$invocation\n";

# Output some more useful information on current settings
if (defined($forced_break)) {
	print "\nWill now proceed to breakdown '$align' using a forced breakpoint after ";
	print "",(($forced_break == 1) ? "every character, " : "every $forced_break characters, "), "and a minimum block size of $min_block_size.\n\n";
}
else {
	print "\nWill now proceed to breakdown '$align' using a minimum block size of $min_block_size.\n\n";
}
print "PAUP settings: gaps will ".((!defined($gap_as_char)) ? "not " : "")."be treated ".
      "as characters, missing and ambiguous sites will ".((defined($exclude_gaps)) ? "not " : "")."be included.\n";
print "MDL settings: nletters = $nletters, nbestpart = $nbestpart, ngroupmax = $ngroupmax.\n\n";

# Relative site numbers of PI characters
my @locations;

# Absolute site numbers of PI characters
my @translated_locations;

my %align = parse_input();
get_informative_chars({'ALIGN' => \%align});
undef(%align);

run_mdl();

sub parse_input {
	my $pipe = new IO::Pipe;
	my $select = new IO::Select; 
	my $pid = fork();

	# By handling the file reading within a fork we save some RAM

	my %align;
	if ($pid == 0) {
		my $TO_PARENT = $pipe->writer();

		# Read first line of input file to determine formatting
		my $first_line;
		open(my $align_file, '<', $align_root) 
			or die "Could not open '$align_root': $!\n";
		while (my $line = <$align_file>) {
			$first_line = $line;
			last;
		}
		close($align_file);

		# Determine whether file is in Nexus or FASTA format, then parse it
		if ($first_line =~ /#NEXUS/i) {
			#print "Input file \"$alignment_path\" appears to be a  Nexus file.\n\n";
			print "Input file '$align' appears to be a  Nexus file.\n";
			%align = parse_nexus($align_root);
		}
		elsif ($first_line =~ /^>/) {
			print "Input file '$align' appears to be a FASTA file.\n";
			%align = parse_fasta($align_root);
		}
		else {
			die "File format for '$align_root' is unrecognized.\n";
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
	my $file_name = shift;

	my $in_data_block;
	my $in_data_matrix;

	my %align;
	open(my $file, "<", $file_name);
	while (my $line = <$file>) {

		# Start of data block
		if ($line =~ /begin data;/i) {
			$in_data_block++;
		}

		# Start of data matrix
		if ($line =~ /matrix/i && $in_data_block) {
			$in_data_matrix++;
		}

		# In data matrix
		if ($in_data_matrix) {
			if ($line =~ /(\S+)\s+(\S+)/) {
				my ($taxon, $sequence) = ($1, $2);

				# Concat allows for reading in multi-line nexus
				$align{$taxon} .= $sequence;
			}
		} 
		
		# End of data matrix
		if ($line =~ /;/ && $in_data_matrix) {
			undef($in_data_matrix);
		}

		# End of data matrix
		if ($line =~ /end;/i && $in_data_block) {
			undef($in_data_block);
		}
	}
	close($file);

	return %align;
}

# Parses alignment data from a properly formatted FASTA file
# Parameters:
#     DATA  - an array reference containing the entire contents of the FASTA file
#     ALIGN - a reference to a hash which will have the alignment data stored with 
#             taxa names as its keys and the nucleotide data as values
sub parse_fasta {
	my $filename = shift;

	my $taxon;
	my %align;
	open(my $alignment_file, '<', $filename) 
		or die "Could not open '$filename': $!\n";

	while (my $line = <$alignment_file>) {
		$line =~ s/^\s+|\s+$//g;

		# Taxon name
		if ($line =~ /^>(.*)/) {
			$taxon = $1;
		}
		else {
			# Taxon sequence
			$taxon =~ s/-/_/g;
			$align{$taxon} .= $line;
		}
	}
	close($alignment_file);
	
	return %align;
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

	# Output data matrix
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

	# File name of log which contains PI character status
	my $log_path = "log-$align_root_no_ext"."-$token-".($number+1)."of$total.nex";

	# Commands we will issue to PAUP
	my $paup_commands = "set monitor=no;\nlog file=$log_path;\nexecute $file_name;\ncstatus full=yes;\nlog stop\n";

	# Temporarily redirect STDOUT to prevent clutter
	open(my $std_out, ">&", *STDOUT);
	close(STDOUT);

	# Run PAUP via a pipe to avoid having to create a command file
	#open(my $paup_pipe, "|-",  $paup, "-L", $log_path, "-n"); 
	open(my $paup_pipe, "|-",  $paup, "-n"); 
	foreach my $command (split("\n", $paup_commands)) {
		print {$paup_pipe} $command,"\n";
	}
	close($paup_pipe);

	# Redirect STDOUT back to terminal
	open(STDOUT, ">&", $std_out);
	close($std_out);

	# Read the resulting log with sysread (faster) 
#	my $log_contents;
#	open(my $log, "<", $log_path) 
#		or die "Could not open '$log_path': $!.\n";
#	sysread($log, $log_contents, -s $log_path);
#	close($log) and unlink($log_path);
#
#	my @log_lines = split("\n", $log_contents);
#
#	# Use map to create an array of only parsimony-informative sites
#	my @filtered_log = grep(/^\d+\s+|(Data matrix)/, @log_lines);
#	my @locations = map { my @line = split(" ", $_);  ($line[2] eq '-' && $line[3] ne 'X') ? $line[0] : () } @filtered_log;
#
#	# Determine total number of characters in the alignment
#	my $total_chars;
#	my $info_line = shift(@filtered_log);
#	if ($info_line =~ /data matrix has \d+ taxa, (\d+) characters/i) { 
#		$total_chars = $1;
#	}
#
#	# Send the locations of the parsimony-informative characters to the parent.
#	foreach my $location (@locations) {
#		print {$$TO_PARENT} "$number-$total_chars--$location\n";
#	}

	my $total_chars;
	open(my $log, "<", $log_path) || die "Could not open '$log_path': $!\n";
	while (my $line = <$log>) {

		# Determine total number of characters in the alignment
		if ($line =~ /data matrix has \d+ taxa, (\d+) characters/i) { 
			$total_chars = $1;
		}
		# Line containing character states
		elsif ($line =~ /^\d+\s+/) {
			my @line = split(" ", $line);	

			# Check that PAUP output character states correctly
			die "Undefined character states detected in PAUP output.\n" if (!defined($line[2]) || !defined($line[3]));

			# Indicates that character is parsimony-informative
			if ($line[2] eq '-' && $line[3] ne 'X') {
				print {$$TO_PARENT} "$number-$total_chars--$line[0]\n";
			}
		}
	}
	close($log);

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

		# Determine number of PICs in each forced breakpoint
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
				# Send this script to the machine
				system("scp", "-q", $script_path, $machine.":/tmp");

				# Send PAUP executable to the machine
				system("scp", "-q", $paup , $machine.":/tmp");

				# Send reduced alignments to remote machines
				system("scp -q $align_root_no_ext-reduced-*.nex $machine:/tmp");

				# Execute this perl script on the given machine
				# -tt forces pseudo-terminal allocation and lets us stop remote processes
				exec("ssh", "-tt", $machine, "perl", "/tmp/$script_name", "--server-ip=$server_ip:$port");
			}
			else {
				# Send this script to the machine
				system("cp", $script_path, "/tmp");

				# Send PAUP executable to the machine
				system("cp", $paup , "/tmp");

				# Execute this perl script on the given machine
				exec("perl", "/tmp/$script_name", "--server-ip=127.0.0.1:$port");
			}

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

		# Free up CPU by sleeping for 10 ms
		usleep(10000);

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
							print "\n  Determining commands for each block... ";
							
							my $time = time();
							write_mdl_command({'GENES' => \@gene_subset});

							print "done. ($start - ".($start + scalar(@gene_subset) - 1).")\n";
							#print "  Total execution time: ", secs_to_readable({'TIME' => time() - $time}), "\n";

							$subset_count++;
						}

						
						# Retrieve the current job and information from the queue
						my $job = shift(@gene_subset);
						my $paup_command = $job->{'PAUP_COMMAND'};
						my $score_file_name = $job->{'SCORE_FILE_NAME'};
						my $align_file_name = $job->{'ALIGN_FILE_NAME'};

						# Check whether the client is remote or local, send it needed files if remote
						if ($client_ip eq $server_ip) {
							print {$client} "CHDIR: ".abs_path($score_dir)."\n";
						}
						print {$client} "NEW: '$score_file_name' '$paup_command'\n";
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
	print "Total execution time: ", sec2human(time() - $time), ".\n\n";

	chdir($project_name);
	%align = parse_input();

	# Calculate number of taxa in alignment for MDL
	my $ntax = scalar(keys %align);

	# Run each partition defined by the forced break points separately through MDL
	foreach my $partition (1 .. $total) {
		my $data_file = "$align_root_no_ext-reduced-$partition"."of$total.nex";
		my $output_name = $partition_dir."$align_root_no_ext-$partition"."of$total.partitions";
		system("$mdl -ntax $ntax -nchar $nchar{$partition} -scorefile 'mdl-scores/$align_root_no_ext-all-scores-$partition' -nletters $nletters -datafile '$data_file' -nbestpart $nbestpart -ngroupmax $ngroupmax -o '$output_name' -ncharbase $min_block_size >/dev/null");
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
	#$command .= "pset gapmode=newstate;\\n" if (!defined($gap_isnt_char));
	$command .= "pset gapmode=newstate;\\n" if (defined($gap_as_char));
	$command .= "include $start-$end / only;\\n";
	$command .= "exclude missambig;\\n" if ($exclude_gaps);
	$command .= "hsearch collapse=no;\\n";
	$command .= "Pscores 1 / scorefile=$score_file_name replace=yes;\\n";
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

	# Recursively run this method until we have the requested number of job descriptions
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

		# Determine each partition's start and end
		# These indices are in parsimony-informative characters
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

	# Check that if we have more than one partition
	if (scalar(keys %reduced_partitions) > 1) {
		foreach my $partition (sort {$a <=> $b} keys %reduced_partitions) {
			my $reduced_start = $reduced_partitions{$partition}->{'START'};
			my $reduced_end = $reduced_partitions{$partition}->{'END'};

			# The first and last partition must be handled slightly differently
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

			# Store breakpoint indices for later output
			$stats{$partition}->{'REDUCED_START'} = $reduced_start;
			$stats{$partition}->{'REDUCED_END'} = $reduced_end;
			$stats{$partition}->{'FULL_START'} = $full_partitions{$partition}->{'START'};
			$stats{$partition}->{'FULL_END'} = $full_partitions{$partition}->{'END'};
		}
	}
	else {
		my $partition = (keys %reduced_partitions)[0];
		my $reduced_start = $reduced_partitions{$partition}->{'START'};
		my $reduced_end = $reduced_partitions{$partition}->{'END'};

		$full_partitions{$partition}->{'START'} = 1;
		$full_partitions{$partition}->{'END'} = length((values %align)[0]); # + 1?

		$stats{$partition}->{'REDUCED_START'} = $reduced_start;
		$stats{$partition}->{'REDUCED_END'} = $reduced_end;
		$stats{$partition}->{'FULL_START'} = $full_partitions{$partition}->{'START'};
		$stats{$partition}->{'FULL_END'} = $full_partitions{$partition}->{'END'};
	}

	print "Alignment has been broken down into ",(scalar(keys %full_partitions))," total partition(s).\n";

	# Output potentially useful information on partitioning
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

	my @reduced_files = glob("$align_root_no_ext-reduced-*.nex");
	print "\nCompressing and archiving parsimony-informative character only alignments... ";
	#system("tar czf $align_root_no_ext-reduced.tar.gz @reduced_files --remove-files");
	system("tar", "czf", "$align_root_no_ext-reduced.tar.gz", @reduced_files);
	unlink(@reduced_files);
	print "done.\n";

	chdir($gene_dir);	

	# Compress and archive results
	my @gene_file_names = glob("$align_root_no_ext-*-*.nex");
	@gene_file_names = sort { local $a = $a; local $b = $b; $a =~ s/.*-(\d+).nex$/$1/; $b =~ s/.*-(\d+).nex$/$1/; $a <=> $b } @gene_file_names;

	# Tarball and zip the resulting sequence partitions
	print "Compressing and archiving partition sequence files... ";
	#system("tar czf $align_root_no_ext.tar.gz @gene_file_names --remove-files");
	system("tar", "czf", "$align_root_no_ext.tar.gz", @gene_file_names);
	unlink(@gene_file_names);
	system("mv '$align_root_no_ext.tar.gz' ..");
	print "done.\n";

	# Tarball and zip the mdl partitioning files
	chdir("..");
	chdir($partition_dir);
	my @partitioning_files = glob("$align_root_no_ext*.partitions*");
	print "Compressing and archiving resulting mdl partitioning output... ";
	#system("tar czf $align_root_no_ext-partitions.tar.gz @partitioning_files --remove-files");
	system("tar", "czf", "$align_root_no_ext-partitions.tar.gz", @partitioning_files);
	unlink(@partitioning_files);
	system("mv '$align_root_no_ext-partitions.tar.gz' ..");
	print "done.\n";

	# Tarball and zip the parismony scores for each possible partition
	chdir("..");
	chdir($score_dir);
	my @score_files = glob("$align_root_no_ext-all-scores-*");
	print "Compressing and archiving individual parsimony scores for each block... ";
	#system("tar czf $align_root_no_ext-scores.tar.gz @score_files --remove-files");
	system("tar", "czf", "$align_root_no_ext-scores.tar.gz", @score_files);
	unlink(@score_files);
	system("mv '$align_root_no_ext-scores.tar.gz' ..");
	print "done.\n";

	# Remove the now empty directories
	chdir("..");
	rmdir($gene_dir);
	rmdir($score_dir);
	rmdir($partition_dir);

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

	# Output file name
	my $file_name = $gene_dir.$align_root_no_ext."-$start-$string_end.nex";
	unlink($file_name) if (-e $file_name);

	sysopen(my $out, $file_name, O_WRONLY|O_NDELAY|O_CREAT)
		or die "Could not open '$file_name': $!.\n";

	my $line_length = 90000;
	my @taxa = sort { $a cmp $b } (keys %{$align});

	print {$out} "#NEXUS\nbegin data;\n dimensions ntax=".scalar(@taxa)." nchar=$nchar;\n ";
	print {$out} "format datatype=dna interleave=yes gap=- missing=?;\n matrix\n";

	# Output the data matrix
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

	# Sanity check to see that each alignment has the proper number of characters
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
	my ($opt_name, $address) = @_;	

	my ($server_ip, $port) = split(":", $address);

	chdir("/tmp");
	my $paup = "/tmp/paup";

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

	# Stores filenames of unneeded files
	my @unlink;

	# Change signal handling so killing the server kills these processes and cleans up
	$SIG{CHLD} = 'IGNORE';
	$SIG{HUP}  = sub { unlink($0, $paup); kill -15, $$; };
	$SIG{TERM} = sub { unlink(@unlink); unlink(glob("*-reduced-*.nex")); exit(0); };

	# Connect to the server
	my $sock = new IO::Socket::INET(
		PeerAddr  => $server_ip.":".$port,
		Proto => 'tcp') 
	or exit(0); 
	$sock->autoflush(1);

	# Ask server for a job
	print {$sock} "NEW: $ip\n";
	while (chomp(my $response = <$sock>)) {

		# Responses to potential server commands
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
			$paup_command =~ s/\.\.\///s if ($server_ip ne "127.0.0.1" && $server_ip ne $ip);

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
	my %actions = ( 'write_partition'           => \&write_partition,
					'write_nexus_file'          => \&write_nexus_file,
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
			my $pid;
			until (defined($pid)) { $pid = fork(); usleep(30000); }
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
			my $pid;
			until (defined($pid)) { $pid = fork(); usleep(30000); }

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

	# Number we want to get the digits of useful for nicely formatted printf output
	my $number = $settings->{'NUMBER'};

	my $digits = 1;
	while (floor($number / 10) != 0) {
		$number = floor($number / 10);
		$digits++;
	}

	return $digits;	
}

sub get_free_cpus {

	# Set to user-specified value if present
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

sub INT_handler {
	clean_up({'DIRS' => 1});
	exit(0);
}

sub clean_up {
	my $settings = shift;

	my $remove_dirs = $settings->{'DIRS'};
	my $current_dir = getcwd();

	chdir($project_name);
	unlink(glob($score_dir."$align_root_no_ext*"));

	if ($remove_dirs) {
		rmdir($gene_dir);
		rmdir($score_dir);
	}
	chdir($current_dir);
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
	return "Usage: mdl.pl [ALIGNMENT FILE] [-b MINIMUM BLOCK LENGTH]\n";
}

sub help {
print <<EOF; 
@{[usage()]}
Use Minimum Description Length (Ané 2011) to break a given alignment into blocks with homologous topologies.

  -b, --block-size       sets the minimum number of parsimony-informative characters which can be found in a block (REQUIRED)
  -f, --forced-break     forces a breakpoint after the specified number of parsimony-informative characters, these partitions 
                         are then run independently this setting is recommended for very large alignments (default: none);
  --gap-as-char          specify to treat gaps as informative characters in PAUP* parsimony analyses
  -o, --out-dir          name of the directory to store output files in (default: "mdl-" + Unix time of script invocation)
  -T, --n-threads        the number of forks ALL hosts running analyses can use concurrently (default: current number of free CPUs)
  --machine-file         file name containing hosts to ssh onto and perform analyses on, passwordless login MUST be enabled
                         for each host specified in this file
  --port                 specifies the port to utilize on the server (Default: 10001)
  --usage                display proper script invocation format
  -h, --help             display this help and exit

Examples:
  perl mdl.pl gene.fa -b 10 --machine-file hosts.txt     runs MDL using computers specified in hosts.txt on gene.fa, breaks can be 
                                                         placed every 10 parsimony-informative characters 
  perl mdl.pl chromosome.fa -b 10 -f 1000                runs MDL on chromosome.fa, a breakpoint will be placed every 1000 parsimony 
                                                         informative characters, within these, further breakpoints can be placed 
                                                         every 10 parsimony-informative characters

Mail bug reports and suggestions to <noah.stenz.github\@gmail.com>
EOF
exit(0);
}
