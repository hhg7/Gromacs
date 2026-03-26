require 5.036;
use 5.042.1;
no source::encoding;
# automatically load to reduce boilerplate
use DDP { output => 'STDOUT', array_max => 10, show_memsize => 1 }; 
use Devel::Confess 'color';
use autodie ':default';
package Gromacs; # packages after this are used, but aren't available when the module is loaded
our $VERSION = 0.01;
use warnings FATAL => 'all';
use autodie ':default';
use File::Temp;
use List::MoreUtils qw(first_index lastidx);
use DDP { output => 'STDOUT', array_max => 10, show_memsize => 1 };
use Devel::Confess 'color';
use Cwd 'getcwd';
use Term::ANSIColor;
use Capture::Tiny 'capture';
use Exporter qw(import);
use List::Util 'first';
use Scalar::Util 'looks_like_number';
use Matplotlib::Simple;
use List::MoreUtils 'firstidx';
use SimpleFlow;
use JSON qw(decode_json encode_json);

our @EXPORT = qw(
add_chain_id2pdb
dssp
gromacs_log2hash
make_rmsd_json
get_group_indices
group_interaction_energy
pdb_coord
pdb_info
plot_gromacs_energy
split_groups
);
our @EXPORT_OK = @EXPORT;

sub json_file_to_ref ($json_filename) {
	die "$json_filename doesn't exist or isn't a file" unless -f $json_filename;
	die "$json_filename has 0 size" if (-s $json_filename == 0);
#	say "Reading $json_filename" if defined $_[0];
	open my $fh, '<:raw', $json_filename; # Read it unmangled
	local $/;                     # Read whole file
	my $json = <$fh>;             # This is UTF-8
#	$json =~ s/NaN/"NaN"/g;
	return decode_json($json); # This produces decoded text
}
sub execute ($cmd, $return = 'exit', $die = 1) {
	if ($return !~ m/^(exit|stdout|stderr|all)$/) {
		die "you gave \$return = \"$return\", while this subroutine only accepts ^(exit|stdout|stderr)\$";
	}
	my ($stdout, $stderr, $exit) = capture {
		system( $cmd )
	};
	if (($die == 1) && ($exit != 0)) {
		say STDERR "exit = $exit";
		say STDERR "STDOUT = $stdout";
		say STDERR "STDERR = $stderr";
		die "$cmd\n failed";
	}
	if ($return eq 'exit') {
		return $exit
	} elsif ($return eq 'stderr') {
		chomp $stderr;
		return $stderr
	} elsif ($return eq 'stdout') {
		chomp $stdout;
		return $stdout
	} elsif ($return eq 'all') {
		chomp $stdout;
		chomp $stderr;
		return {
			exit   => $exit, 
			stdout => $stdout, 
			stderr => $stderr
		}
	} else {
		die "$return broke pigeonholes"
	}
	return $stdout
}

sub ref_to_json_file ($ref, $json_filename) {
	my $ref_json_filename = ref $json_filename;
	unless ($ref_json_filename eq '') {
		die "$json_filename isn't a scalar/string";
	}
	open my $fh, '>:raw', $json_filename; # Write it unmangled
	say $fh encode_json($ref);
	say 'Wrote ' . colored(['blue on_red'], $json_filename);
	return $json_filename;
}

sub gromacs_log2hash ($log, $json_filename = undef) {
	open my $fh, '<', $log;
	my @log = <$fh>;
	close $fh;
	my (@energies, %data, @splice);
	foreach my $i ( grep {$log[$_] =~ m/^Finished\h/} reverse ( ($#log - 9) .. $#log ) ) {
		chomp $log[$i];
		$data{'Finished'} = $log[$i];
		splice @log, $i, $#log - $i; # remove this line and everything after it
		last;
	}
	foreach my $i (reverse 0..$#log) { # remove blank lines at the end
		last if $log[$i] =~ m/\S/;
		pop @log;
	}
	chomp @log;
	my $i = lastidx { /^Performance:\h+\d+/ } @log; # searching from the end to be faster
	if (
			($i > 0)
			&&
			($log[$i] =~ m/^Performance:\h+(\d+)\.(\d+)\h+(\d+)\.(\d+)/)
		) {
		my @metrics = ('ns/day', 'hour/ns', 'ms/step', 'Matom*steps/s');
		my @vars = split /\h+/, $log[$i];
		shift @vars; # skip "Performance:"
		foreach my $metric (@metrics) {
			$data{Performance}{$metric} = shift @vars;
			last if scalar @vars == 0;
		}
		splice @log, $i-1, 2;
	}
	foreach my $i (reverse 0..$#log) {
		if ($log[$i] =~ m/^\h+Time:\h+(\d+)\.(\d+)\h+(\d+)\.(\d+)\h+(\d+)\.(\d+)/) {
			$data{Performance}{'Core Time (s)'} = "$1.$2";
			$data{Performance}{'Wall Time (s)'} = "$3.$4";
			$data{Performance}{'Time (%)'} = "$5.$6";
			splice @log, $i, 1;
			$data{Performance}{'Wall Time (min)'} = $data{Performance}{'Wall Time (s)'} / 60;
			$data{Performance}{'Wall Time (hr)'}  = $data{Performance}{'Wall Time (s)'} / 3600;
			last;
		}
	}
	$i = lastidx {/^\h*Breakdown of PME mesh activities/} @log;
	if ($i > 0) {
		my @j = grep {$log[$_] =~ m/^\-+$/} $i..$#log;
		die "no lines with \"---...\" found in $log" if scalar @j == 0;
		my @col = (
		#'Activity', # activity is the hash key, which isn't used in the hash slice
		'Num Ranks', 'Num Threads', 'Call Count', 'Wall time (s)', 'Giga-Cycles total', 'Giga-Cycles %');
		foreach my $line ($j[0]+1 .. ($j[1] - 1)) {
			$log[$line] =~ s/^\h+//;
			chomp $log[$line];
			my @line = grep {$_ ne ''} split /\h{2,}/, $log[$line];
			my $key = shift @line;
			if (scalar @line == scalar @col) {
				@{ $data{'REAL CYCLE AND TIME ACCOUNTING'}{$key} }{@col} = @line;
			} else {
				p $log[$line];
				die "The above line in $log failed pigeonholes.";
			}
		}
		splice @log, $i, $i+$j[0];
	}
	$i = lastidx {/^\h*Activity:\h/} @log;
	if ($i > 0) {
		my @j = grep {$log[$_] =~ m/^\-+$/} $i..$#log;
		die "no lines with \"---...\" found in $log" if scalar @j == 0;
		my @col = (
		#'Activity', # activity is the hash key, which isn't used in the hash slice
		'Num Ranks', 'Num Threads', 'Call Count', 'Wall time (s)', 'Giga-Cycles total', 'Giga-Cycles %');
		my @wait_col = @col[3..5];
		foreach my $line ($j[0]+1 .. ($j[1] - 1)) {
			$log[$line] =~ s/^\h+//;
			chomp $log[$line];
			my @line = grep {$_ ne ''} split /\h{2,}/, $log[$line];
			my $key = shift @line;
			if (scalar @line == scalar @col) {
				@{ $data{'REAL CYCLE AND TIME ACCOUNTING'}{$key} }{@col} = @line;
			} elsif (scalar @line == scalar @wait_col) {
				@{ $data{'REAL CYCLE AND TIME ACCOUNTING'}{$key} }{@wait_col} = @line;
			} else {
				p $log[$line];
				die "The above line in $log failed pigeonholes.";
			}
		}
		splice @log, $i, $i+$j[0];
	}
	$i = lastidx {/^\h*Computing:\h/} @log;
	if ($i > 0) { # -1 if it doesn't show up with "first_index"
		my @j = grep {$log[$_] =~ m/^\-+$/} $i..$#log;
		die "no lines with \"---...\" found in $log" if scalar @j == 0;
		$log[$i] =~ s/^\h+//;
		my @col = split /\h{2,}/, $log[$i];
		shift @col;
		foreach my $line ($j[0]+1 .. ($j[1]-1)) {
			my $original_line = $log[$line];
			$log[$line] =~ s/^\h+//;
			my @line = grep {$_ ne ''} split /\h{2,}/, $log[$line];
			my $key = shift @line;
			$log[$line] =~ s/^\Q$key\E\h*//;
			@line = split /\h+/, $log[$line];
			if (scalar @line != (scalar @col)) {
				say STDERR 'For these columns:';
				p @col;
				say STDERR 'This line is no good:';
				p @line;
				die "these columns do not line up for \"$original_line\" in $log";
			}
			$key =~ s/^\h+//;
			@{ $data{'MEGA-FLOPS ACCOUNTING'}{$key} }{@col} = @line; # hash slice
		}
		my @line = grep {$_ ne ''} split /\h+/, $log[$j[1]+1];
		if ($line[0] ne 'Total') {
			die "\"Total\" not found for Mega-flops accounting in $log";
		}
		@{ $data{'MEGA-FLOPS ACCOUNTING'}{Total} }{'M-Flops', '% Flops'} = @line[1,2];
		splice @log, $i, $j[1] - $i;
	}
	foreach my $i (grep {$log[$_] =~ m/^DD\h+step\h+\d+/} reverse 0..$#log) {
		splice @log, $i-1, 2; # remove this line and the blank that comes before it
	}
	foreach my $i (grep {$log[$_] =~ m/^step\h+\d+\h+Turning on/} reverse 0..$#log) {
		splice @log, $i, 1;
	}
	foreach my $i (grep {$log[$_] =~ m/^Current ref_t for group System:\h+\d+\.\d+/} reverse 0..$#log) {
		splice @log, $i-1, 2;
	}
	foreach my $i (grep {$log[$_] =~ m/^\++\hPLEASE READ AND CITE THE FOLLOWING REFERENCE\h\++$/} 0..$#log) {
		my $thanks_line;
		foreach my $j ($i..$#log) {
			if ($log[$j] =~ m/^[\-\h]+.+Thank You\h[\h\-]+$/) {
				$thanks_line = $j;
				last;
			}
		}
		if (not defined $thanks_line) {
			die "Couldn't find \"Thank You\" line in $log after \"$log[$i]\"";
		}
		push @splice, [$i, $thanks_line - $i];
	}
	foreach my $spl (reverse @splice) {
		splice @log, $spl->[0], $spl->[1];
	}
	$i = first_index {/^qm-opts:$/} @log;
	if ($i > 0) {
		my @j = grep {$log[$_] =~ /^\H/} $i+1..$#log;
		foreach my $l ($i+1..($j[0]-1)) {
			next if $log[$l] !~ m/^   \w/;
			next if $log[$l] !~ m/\h=\h/;
			$log[$l] =~ s/^\h+//;
			my @line = split /\h+=\h+/, $log[$l];
			$data{'qm-opts'}{$line[0]} = $line[1];
		}
		# remove these lines to save time and prevent errors
		splice @log, $i, $j[0] - $i;
	}
	my @colvars_lines = grep {/^colvars:\h/} @log;
	if (scalar @colvars_lines > 0) {
		my $last_index = first_index {/^colvars:\h+The final output state file will be ".+"/} @colvars_lines;
		if ($last_index > 0) {
			splice @colvars_lines, $last_index, scalar @colvars_lines - 1 - $last_index;
		}
		my @colvars_start_idx = grep {$colvars_lines[$_] eq 'colvars:   Initializing a new collective variable.'} 0..$#colvars_lines;
		foreach my ($idx, $start_idx) (indexed @colvars_start_idx) {
			my $end_idx = first {
									$colvars_lines[$_] =~ m/^colvars:\h+\-+$/
								} $start_idx+1..$#colvars_lines;
			unless (defined $end_idx) {
				die "Couldn't find an end for the colvar after index $start_idx in $log";
			}
			my @colvar = @colvars_lines[$start_idx+1..$end_idx];
			if (scalar @colvar == 0) {
				die "0 elements in array: Couldn't get colvar data from $log";
			}
			my $name;
			foreach my ($idx, $line) (indexed @colvar) {
				if ($line =~ m/^colvars:\h+#\h+name = "([^"]+)/) {
					$name = $1;
					splice @colvar, $idx, 1;
					next;
				}
			}
			unless (defined $name) {
				p @colvar;
				die "Couldn't get name for colvar in $log";
			}
			# get info for separate atom groups
			my @atom_grp_start = grep {$colvar[$_] =~ m/^colvars:\h+Initializing atom group/} 0..$#colvar;
			if (scalar @atom_grp_start == 0) {
				p @colvar;
				die "Couldn't get atom group starts from the above colvar in $log";
			}
			my @atom_grp_end   = grep {$colvar[$_] =~ m/^colvars:\h+Atom group \".+\" defined/} ($atom_grp_start[0]+1)..$#colvar;
			if (scalar @atom_grp_start != scalar @atom_grp_end) {
				p @colvar;
				die "The above colvar group does not have the same # of starting and ending lines in $log";
			}
			foreach my ($g, $atom_grp_start) (indexed @atom_grp_start) {
				my $grp_name;
				foreach my ($ag, $line) (indexed @colvar[$atom_grp_start..$atom_grp_end[$g]]) {
					if ($line =~ m/^colvars:\h+Initializing atom group "(.+)"/) {
						$grp_name = $1;
						next;
					}
					unless (defined $grp_name) {
						p @colvar;
						die "No group name is defined for group $g in $log";
					}
					$line =~ s/^colvars:\h+#?\h*//;
					if ($line =~ m/Atom group \".+\" defined with (\d+) atoms requested: total mass = (\d+)\.(\d+), total charge = (\-?)([\d+\.]+)/) {
						$data{colvars}{$name}{$grp_name}{'No. of atoms'} = $1;
						$data{colvars}{$name}{$grp_name}{'total mass'}   = "$2.$3";
						$data{colvars}{$name}{$grp_name}{'total charge'} = "$4$5";
						$data{colvars}{$name}{$grp_name}{'total charge'} =~ s/\.$//;
						last;
					}
					my @line = split /\h+=\h+/, $line, 2;
					if (scalar @line != 2) {
						p @line;
						die "$log has a line with no. of elements != 2; where \$ag = $ag within [$atom_grp_start..$atom_grp_end[$g]]";
					}
					if ($line[1] =~ m/^"(.+)"$/) {
						$line[1] = $1;
					}
					$data{colvars}{$name}{$grp_name}{$line[0]} = $line[1];
				}
			}
			@atom_grp_start = reverse @atom_grp_start;
			@atom_grp_end   = reverse @atom_grp_end;
			foreach my ($grp, $start) (indexed @atom_grp_start) {
				# remove these lines so they're not double-counted
				splice @colvar, $start, $atom_grp_end[$grp] - $start; 
			}
			foreach my $line (grep {/^colvars:\h+#\h/} @colvar) {
				$line =~ s/^colvars:\h+#\h+//;
				if ($line =~ m/^Initializing a new "(.+)" component/) {
					$data{colvars}{$name}{type} = $1;
					next;
				}
				my @line = split /\h+=\h+/, $line;
				if (scalar @line != 2) {
					p @line;
					die "$log has a line with no. of elements != 2";
				}
				foreach my $i (@line) {
					$i =~ s/^\h+//;
					$i =~ s/\h+$//;
				}
				$data{colvars}{$name}{$line[0]} = $line[1];
			}
		}
# ensure that I got the correct number of collective variables
		my $no_gromacs_colvars;
		foreach my $line (grep {
						/^colvars:\h+Collective variables initialized, (\d+) in total\.$/
									} @colvars_lines) {
			if ($line =~ m/(\d+) in total\.$/) {
				$no_gromacs_colvars = $1;
			} else {
				die "$line failed regex.";
			}
			last;
		}
		my $no_detected_colvars = scalar keys %{ $data{colvars} };
		if ($no_gromacs_colvars != $no_detected_colvars) {
			p $data{colvars};
			die "Gromacs wrote that there are $no_gromacs_colvars colvars, but I gathered $no_detected_colvars from $log";
		}
		@colvars_start_idx = grep {$colvars_lines[$_] =~ m/^colvars:\h+Initializing a new ".+" instance\.$/} 0..$#colvars_lines;
		foreach my $start_idx (@colvars_start_idx) {
			my $instance;
			if ($colvars_lines[$start_idx] =~ m/^colvars:\h+Initializing a new "(.+)" instance\.$/) {
				$instance = $1;
			} else {
				die "$colvars_lines[$start_idx] from $log failed regex";
			}
			my $end_idx = first {
									$colvars_lines[$_] =~ m/^colvars:\h+\-+$/
								} $start_idx..$#colvars_lines;
			unless (defined $end_idx) {
				die "Couldn't find an end for the colvar after index $start_idx in $log";
			}
			foreach my $line ( @colvars_lines[$start_idx+1..$end_idx-1] ) {
				$line =~ s/^colvars:\h+#\h*//;
				my @line = split /\h*=\h*/, $line;
				foreach my $i (@line) {
					$i =~ s/\{//;
					$i =~ s/}//;
					$i =~ s/^\h+//;
					$i =~ s/\h+$//;
				}
				$data{colvars}{$instance}{$line[0]} = $line[1];
			}
		}
		my $i = firstidx {/^colvars:\h+The following index groups are currently defined:/} @colvars_lines;
		my $last_i = first {$colvars_lines[$_] =~ m/^colvars:\h+#/} ($i+1)..$#colvars_lines;
		foreach my $i (($i+1)..($last_i-1)) {
			if ($colvars_lines[$i] =~ m/^colvars:\h+(.+) \((\d+)/) {
				$data{colvars}{'index group # of atoms'}{$1} = $2;
			} else {
				die "$colvars_lines[$i]: Failed regex in $log";
			}
		}
	}
	@log = grep {$_ !~ m/^colvars:\h/} @log; # I've already dealt with these, remove to prevent bugs
	$i = first_index {/^\h+CPU info:/} @log;
	if ($i > 0) {
		my @j = grep {$log[$_] eq ''} ($i+1)..$#log;
		foreach my $j ($i+1..$j[0]-1) {
			next if $log[$j] !~ m/:.+/; # ":" cannot end the line
			$log[$j] =~ s/^\h+//;
			my @line = split /\h*:\h*/, $log[$j];
			if (scalar @line == 2) {
				$data{'CPU info'}{$line[0]} = $line[1];
			}
		}
		splice @log, $i, $j[0] - $i;
	}
	foreach my $line (grep {/^There are: \d+ Atoms$/} @log) {
		if ($line =~ m/(\d+)/) {
			$data{'Atom Count'} = $1;
			last;
		} else {
			die "$line failed regex.";
		}
	}
	foreach my $writing_index (reverse grep { $log[$_] =~ m/^Writing checkpoint/} 99..$#log) {
		splice @log, $writing_index, 3;
	}
	my $input_param_i = first_index {$_ eq 'Input Parameters:'} @log;
	my $ending_i = first {$log[$_] =~ m/^\S/} $input_param_i+1..$#log;
	foreach my $i (reverse $input_param_i..$ending_i) {
		$log[$i] =~ s/^\h+//;
		my @line = split /\h+=\h+/, $log[$i];
		next unless scalar @line == 2;
		$data{'Input Parameters'}{$line[0]} = $line[1];
	}
	splice @log, $input_param_i, $ending_i - $input_param_i; # prevent regex confusion later on
	splice @log, 0, first_index { /^GROMACS:\h+/ } @log;
	my @time_indices = grep {
									$log[$_-1] =~ m/^\h+Step\h+Time$/
								&&
									$log[$_] =~ m/
									^\h+
									\d+        # Step
									\h+
									[\d\.]+    # Time (10.000)
									$/x
								&&
									$log[$_+1] eq ''
								&&
									(
										($log[$_+2] =~ m/^\h+Energies \(kJ\/mol\)/)
										||
										($log[$_+3] =~ m/^\h+Energies \(kJ\/mol\)/)
									)
									
								} 0..$#log-3;
	if (scalar @time_indices == 0) {
		p @log, array_max => scalar @log;
		my @step_time_i = grep {$log[$_] =~ m/^\h+Step\h+Time$/} 0..$#log;
		my @t = @log[@step_time_i];
		p @t;
		die "Couldn't get times for $log";
	}
	foreach my $time_index (@time_indices) {
		if ($log[$time_index] =~ m/(\d+)\.(\d+)$/) {
			push @{ $data{time} }, "$1.$2";
		} else {
			die "$log[$time_index] failed regex.";
		}
	}
	my $reading_energies = 0;
	foreach (@log) {
		if ($_ eq '   Energies (kJ/mol)') {
			$reading_energies = 1;
			next;
		}
		if ($_ eq '') {
			$reading_energies = 0;
			next;
		}
		next unless $reading_energies == 1;
		if (/^\h+[A-Z]/) {
			while ($_ =~ m/(.{1,15})/g) {
				my $e = $1;
				$e =~ s/^\h+//;
				push @energies, $e;
			}
			next;
		}
		if ((/^\h+\-?\d/) && (scalar @energies > 0)) {
			$_ =~ s/^\h+//;
			my @line = split /\s+/, $_;
			if (scalar @line != scalar @energies) {
				p @energies;
				p @line;
				die "$log line $. has " . scalar @line . ' energies, but should have ' . scalar @energies . "\n";
			}
			foreach my $energy (@energies) {
				push @{ $data{Energy}{$energy} }, shift @line;
			}
			undef @energies;
		}
	}
	if (scalar keys %{ $data{Energy} } == 0) {
		die "Couldn't get any energies from $log";
	}
	foreach my $line (grep {$_ =~ m/^.+:\h+\H+/} @log) {
		my @line = split /:\h+/, $line;
		$data{$line[0]} = $line[1];
	}
	foreach my $key (grep {/^\h+/} keys %data) {
		my $new_key = $key;
		$new_key =~ s/^\h+//;
		$data{$new_key} = delete $data{$key};
	}
	foreach my $i (0..$#log) {
		if ($log[$i] =~ m/^Running on (\d+) nodes? with total (\d+) cores, (\d+) processing units$/) {
			$data{nodes} = $1;
			$data{cores} = $2;
			$data{'processing units'} = $3;
			last;
		}
#		splice @log, $i, 1;
	}
	foreach my $i (0..$#log) {
		if ($log[$i] =~ m/^Hardware detected on host (.+):$/) {
			$data{Host} = $1;
			last;
		}
#		splice @log, $i, 1; # don't possibly confuse future regex
	}
	return \%data unless defined $json_filename;
	open $fh, '>:raw', $json_filename; # Write it unmangled
	say $fh encode_json(\%data);
	say 'Wrote ' . colored(['blue on_red'], $json_filename);
	return \%data;
}

sub plot_gromacs_energy (
	$gromacs,      # a file or a hash from gromacs_log2hash
	$output_file,  # output file image ("file.svg")
	$suptitle      # sup title for output image
) {
	if (ref $gromacs eq '') { # i.e. this is a file
		$gromacs = gromacs_log2hash( $gromacs );
	}
	my @plots;
	foreach my $energy (sort keys %{ $gromacs->{Energy} }) {
		my $ylab = '(kJ/mol)';
		$ylab = 'bar' if $energy =~ m/^Pres/;
		$ylab = 'K'   if $energy eq 'Temperature';
		$ylab = 'Å'   if $energy =~ m/RMSD/i;
		push @plots, {
			data        => {
				$energy => [
					$gromacs->{'time'},
					$gromacs->{Energy}{$energy}
				]
			},
			'plot.type'   => 'plot',
			'show.legend' => 0,
			title         => $energy,
			ylabel        => $ylab,
		};
	}
	plt({
		'output.file' => $output_file,
		plots         => \@plots,
		ncols         => 3,
		nrows         => scalar @plots / 3,
		scale         => 2,
		suptitle      => $suptitle
	});
}

sub make_rmsd_json ($xtc, $tpr, $stem) {
	my $datfile = "$stem.rmsd.dat";
	task({
#		'dry.run'      => 1,
		'input.files'  => [$xtc, $tpr],
		cmd            => "printf \"0\\n0\\n\" | gmx rms -s $tpr -f $xtc -m $stem.rmsd.xpm -bin $datfile",
		'output.files' => ["$stem.rmsd.xpm", $datfile]
	});
	my $rmsdjson = $datfile;
	$rmsdjson =~ s/dat$/json/;
	my $py = File::Temp->new(DIR => '/tmp', SUFFIX => '.py', UNLINK => 1);
	say $py 'import numpy as np
import json
import sys

def ref_to_json_file(data, filename):
	json1=json.dumps(data)
	f = open(filename,"w+")
	print(json1,file=f)';
# 1. Read the binary data file
# The data is a raw dump of 32-bit floats (np.float32)
	say $py "binary_file = '$datfile' # Assuming you named the output with -bin as .dat";
	say $py 'data = np.fromfile(binary_file, dtype=np.float32)';
	say $py 'num_frames = int(np.sqrt(data.size))

if num_frames * num_frames != data.size:
	sys.exit("Error: Data size does not correspond to a square matrix. Check your input.")
else:# 3. Reshape the 1D array into the N x N matrix
	rmsd_matrix = data.reshape(num_frames, num_frames)';
	say $py	"ref_to_json_file(rmsd_matrix.tolist(), '$rmsdjson')";
	close $py;
	task({
		cmd            => 'python ' . $py->filename,
#		'dry.run'      => 1,
		'input.files'  => [$py->filename],
		'output.files' => [$rmsdjson],
	});
}

sub get_group_indices ($ndx_file) {
# purpose is to take the below table, and return an array, where Nr. is the index for "Group"
#Contents of index file cpx.ndx
#--------------------------------------------------
#Nr.   Group               #Entries   First    Last
#   0  System                 56416       1   56416
#   1  Protein                 1136       1    1136
#   2  Receptor                 965       1     965
#   3  Ligand_CA                 10     970    1127
#   4  Ligand                   171     966    1136
#   5  BindingSite_CA            15      58     679
#   6  BindingSite              206      56     700
#   7  Receptor-BindingSite     640       1     965
#   8  BindingAxisAnchor         30     283     702
#----------
# returns a hash
#----------
# BindingAxisAnchor      8,
# BindingSite            6,
# BindingSite_CA         5,
# Ligand                 4,
# Ligand_CA              3,
# Protein                1,
# Receptor               2,
# Receptor-BindingSite   7,
# System                 0

	my $t = task({
		cmd           => "gmx check -n $ndx_file",
		'input.files' => $ndx_file,
		note          => 'read ndx file'
	});
	my @ndx = split /\n/, $t->{stdout};
	my $i = first_index {/^\h+\d+\h/} @ndx;
	splice @ndx, 0, $i;
	my %group2i;
	foreach my $line (@ndx) {
		my @line = grep {$_ =~ m/\H/} split /\h+/, $line, 4;
		$group2i{$line[1]} = $line[0];
	}
	return \%group2i;
}

sub group_interaction_energy ($in_xtc, $out_xtc) {
# this removes the non-Protein elements to get energy between protein and ligand/peptide in MD run
	return 1 if ((-f $out_xtc) && (-s $out_xtc > 9));
	my $progro = '2eq_out.noSOL.gro';
	task({ # group 1 = "Protein"
		cmd            => "echo 1 | gmx trjconv -f 2eq_out.gro -s 2eq_out.gro -o $progro -n cpx.ndx",
		'input.files'  => ['cpx.ndx', '2eq_out.gro'],
		'output.files' => $progro,
		overwrite      => 0
	});
	my $noSOL_top = 'cpx.noSOL.top';
	unless (-f $noSOL_top) {
		open my $in, '<', 'cpx.top';
		open my $out, '>', $noSOL_top;
		while (<$in>) {
			next if /^(?:SOL|NA|CL)\h/;
			next if /^Ion_chain/;
			print $out $_;
		}
	}
	my $mdp = '3md.rerun.mdp';
	unless (-f $mdp) {
		open my $fh, '>', $mdp;
		say $fh 'energygrps = Receptor Ligand';
	}
	my $tpr = '3md.rerun.tpr';
	task({
		cmd            => "gmx grompp -f $mdp -c $progro -p $noSOL_top -n cpx.ndx -o $tpr",
		'input.files'  => [$progro, $noSOL_top, 'cpx.ndx'],
		'output.files' => $tpr,
		overwrite      => 0
	});
	my $stem = $out_xtc;
	$stem =~ s/\.xtc$//;
	task({
		'input.files' => [$tpr, $in_xtc],
		cmd           => "gmx mdrun -s $tpr -rerun $in_xtc -e $stem.edr -nb cpu -deffnm $stem",
		'output.files'=> "$stem.log"
	});
	return $out_xtc;
}
sub dssp ($xtc_files, $pdb_files, $json_files) { # get secondary structure prediction
# Ligand.pdb is made thus: echo 4 | gmx editconf -f cpx.pdb -n cpx.ndx -o Ligand.pdb
	my $n_xtc  = scalar @{ $xtc_files };
	my $n_pdb  = scalar @{ $pdb_files };
	my $n_json = scalar @{ $json_files };
	unless ($n_xtc == $n_pdb == $n_json) {
		die "There are $n_xtc xtc, $n_pdb PDB, and $n_json JSON files.  The list sizes must all be equal";
	}
	my @not_files = grep {not -f -r $_} @{ $xtc_files }, @{ $pdb_files };
	if (scalar @not_files > 0) {
		p @not_files;
		die 'the above files are not files';
	}
	my @c = caller;
	my $py = File::Temp->new(DIR => '/tmp', SUFFIX => '.py', UNLINK => 1);
	say $py 'import mdtraj as md';
	say $py 'import json';
	say $py 'def ref_to_json_file(data, filename):';
	say $py '	json1=json.dumps(data)';
	say $py '	f = open(filename,"w+")';
	say $py '	print(json1,file=f)';
	foreach my ($idx, $xtc) (indexed @{ $xtc_files }) {
		say $py "traj = md.load('$xtc', top = '$pdb_files->[$idx]')";
		say $py 'dssp = md.compute_dssp(traj, simplified = False)';
		say $py 'dssp = dssp.transpose()';
		say $py 'dssp = dssp.tolist()'; # can't export to JSON otherwise
		say $py "ref_to_json_file(dssp, '$json_files->[$idx]')";
	}
	close $py;
	task({
		cmd            => 'python ' . $py->filename,
		'input.files'  => $xtc_files,
		note           => "extracting DSSP from XTC files from $c[1] line $c[2]",
		'output.files' => $json_files,
		overwrite      => 1
	});
}
sub file2string ($file) {
	open my $fh, '<', $file;
	return do { local $/; <$fh> };
}
sub pdb_info (
	$id,        # 2puy
	$stem = $id,# usually the id itself
	$redo = 0   # redo or not
) {
	my $current_sub = (split(/::/,(caller(0))[3]))[-1];
# return a hash of PDB information
# saves to $id.json and $id.html in the directory that the subroutine is called from
	$id = lc $id; # make all characters lower case
	my $json_file = "$stem.json";
	unlink $json_file if (($redo) && (-f $json_file));
	if (-f $json_file) {
		return json_file_to_ref( $json_file );
	}
	my ($html, $json);
	if (-f $json_file) {
		$json = file2string($json_file);
	} else {
		$json = execute("curl -X GET 'https://data.rcsb.org/rest/v1/core/entry/$id'", 'stdout');
		open my $j, '>', $json_file;
		say $j $json;
	}
	$json =~ s/NaN/***NaN***/g;
	$json = decode_json( $json );
	delete $json->{entry}; # worthless, ID is in the name
	delete $json->{rcsb_id} if defined $json->{rcsb_id};
	delete $json->{rcsb_primary_citation}{rcsb_orcididentifiers} if defined $json->{rcsb_primary_citation}{rcsb_orcididentifiers};
	delete $json->{rcsb_entry_container_identifiers}{entry_id} if defined $json->{rcsb_entry_container_identifiers}{entry_id};
	delete $json->{rcsb_entry_container_identifiers}{rcsb_id}  if defined $json->{rcsb_entry_container_identifiers}{rcsb_id};
	my $html_file = "$stem.html";
	if (-f $html_file) {
		$html = file2string( $html_file );
	} else {
		
		$html = execute("wget --quiet --output-document=- https://www.rcsb.org/structure/$id", 'stdout');
		open my $h, '>', $html_file;
		say $h $html;
		say "Wrote $html_file";
	}
	die "$html_file doesn't exist" unless -f $html_file;
	my @species = ($html =~ m/
		rcsb_entity_source_organism\.taxonomy_lineage\.name:([^"]+)
	/xxg);
	@{ $json->{Organisms} } = do { my %seen; grep { !$seen{$_}++} @species }; # get unique species
	foreach my $r ('Resolution', 'R-Value Free', 'R-Value Work', 'R-Value Observed') {
		@{ $json->{$r} } = ($html =~ m/$r:&nbsp<\/strong>([\d\.]+)/g);
		next if scalar @{ $json->{$r} } == 0;
		while ( scalar @{ $json->{$r} } > (scalar @{ $json->{exptl} })) {
			pop @{ $json->{$r} }; # sometimes written 2x on the HTML
		}
	}
	@{ $json->{'Gene Names'} } = ($html =~ m/
		rcsb_entity_source_organism.rcsb_gene_name.value:([^"]+)
	/xxg);
	if ($html =~ m/Mutation\(s\):&nbsp<\/strong>(Yes|No)/) {
		$json->{mutation} = $1;
	} else {
		say "Couldn't find mutation for $id";
	}
	if ($html =~ m/ELECTRON MICROSCOPY<\/li><li id="exp_header_0_em_resolution"><strong>Resolution:&nbsp<\/strong>[\d\.]+&nbsp&Aring;<\/li><li id="exp_header_0_em_aggregationState"><strong>Aggregation State:&nbsp<\/strong>([A-Z\d\s]+)/) {
		$json->{'Aggregation State'} = $1;
	}
=if ($html =~ m/X-RAY DIFFRACTION<\/li><li id="exp_header_0_diffraction_resolution"><strong>Resolution:&nbsp<\/strong>([\d\.]+) Å<\/li><li id="exp_header_0_diffraction_rvalueWork"><strong>R-Value Work:&nbsp<\/strong>([\.\d]+)/) {
		$json->{Resolution}     = $1;
	}
	if (
			(defined $json->{exptl}) &&
			((ref $json->{exptl}) eq 'ARRAY') &&
			((ref $json->{exptl}[0]) eq 'HASH')
		) {
		$json->{Method} = $json->{exptl}[0]{method};
	}
=cut
	if ($html =~ m/>Reconstruction Method:&nbsp<\/strong>([A-Z\s]+)/) {
		$json->{'Reconstruction Method'} = $1;
	}
=if (
			(not defined $json->{Resolution})   &&
			($html =~ /(\d+).(\d+) Å resolution/)
		) {
		$json->{Resolution} = "$1.$2"; # 7o2w, the resolution is written inside the HMTL, but not visible on the webpage
	}
	if (
			(not defined $json->{Resolution})   &&
			#           Resolution:&nbsp</strong>2.10 Å<
			($html =~ m/Resolution:&nbsp<\/strong>([\.\d]+)/)#&nbsp&Aring/)
		) {
		$json->{Resolution} = $1; # 7wzw
	}
=cut
	if ($html =~ m/Aggregation State:&nbsp<\/strong>([A-Z\d]+)/) {
		$json->{'Aggregation State'} = $1;
	}
	if ($html =~ m/<title>([^>]+)<\/title>/) {
		$json->{title} = $1;
	}
	if ($html =~ m/Total Structure Weight: ([\d\.,]+) ([A-Za-z]+)&/) {
		$json->{'Total Structure Weight'} = "$1 $2";
	}
	foreach my $q ('Atom Count', 'Modelled Residue Count', 'Deposited Residue Count') {
		if ($html =~ m/$q: ([\d,]+)/) {
			$json->{$q} = $1;
		}
	}
	my @undef = grep {not defined $json->{$_}} ('Resolution', 'Aggregation State', 'Total Structure Weight', 'Atom Count', 'Modelled Residue Count', 'Deposited Residue Count');
	if (($id ne '7q2y') && ($html =~ m/ELECTRON MICROSCOPY/) && (scalar @undef > 0)) {
		p $html;
		p $json;
		say STDERR "the arguments below for ID = $id are missing";
		p @undef;
		die "Couldn't get a method, resolution, and/or aggregation state for $id, the missing methods are above";
	}
	ref_to_json_file($json, $json_file);
	return $json;
}
sub pdb_coord ($line) {
# https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html#note5
	my %r;
	if ($line =~ m/^(ATOM|HETATM|TER)/) {
		$r{'Record Type'} = $1;
	} else {
		die "$line isn't defined";
	}
	chomp $line;
	my @line = split '', $line;
	if ($r{'Record Type'} =~ m/^(?:ATOM|HETATM)$/) {
		$r{'Atom serial number'}           = join ('', grep {/\d/} @line[6..10]);
		unless (looks_like_number($r{'Atom serial number'})) {
			say STDERR "detected serial number: $r{'Atom serial number'}";
			say STDERR $line;
			die 'Atom serial number is not numeric';
		}
		$r{'Atom name'}                    = join ('', grep {/\H/} @line[12..15]);
		$r{'Alternate location indicator'} = $line[16];
		$r{'Residue name'}                 = join ('', grep {/\H/} @line[17..19]);
		$r{'Chain identifier'}             = $line[21];
		$r{'Residue sequence number'}      = join ('', grep {/\d/} @line[22..25]);
		$r{'Code for insertions of residues'} = $line[26];
		$r{'X orthogonal Å coordinate'}    = join ('', grep {/\H/} @line[30..37]);
		$r{'Y orthogonal Å coordinate'}    = join ('', grep {/\H/} @line[38..45]);
		$r{'Z orthogonal Å coordinate'}    = join ('', grep {/\H/} @line[46..53]);
		$r{'Occupancy'}                    = join ('', grep {/\H/} @line[54..59]);
		$r{'Temperature factor'}           = join ('', grep {/\H/} @line[60..65]);
		$r{'Segment identifier'}           = join ('', grep {/\H/} @line[72..75]);
		$r{'Element symbol'}               = join ('', grep {/\H/} @line[76,77]);
		if ((defined $line[79]) && (grep {/./} @line[78,79])) {
			$r{Charge}                       = join ('', grep {/\H/} @line[78,79]);
		}
	} elsif ($r{'Record Type'} eq 'TER') {
		$r{'Serial Number'} = join ('', grep {/\d/} @line[6..10]);
		$r{'Residue Name'}  = join ('', grep {/\H/} @line[17..19]);
		$r{'Chain identifier'} = $line[21];
		$r{'Residue sequence number'} = join ('', grep {/\d/} @line[22..25]);
		$r{'Code for insertions of residues'} = $line[26];
	}
	return \%r;
}
sub add_chain_id2pdb (
# trjconv omits chain IDs, which makes PyMol visualization look bad
	$pdb_file_with_chains,
	$pdb_file_without_chains,
	$output_file
) {
	my %res_num2chain;
	open my $in, '<', $pdb_file_with_chains;
	while (<$in>) {
		next unless (/^(?:ATOM|HETATM|TER)/);
		my $line = pdb_coord($_);
		my @undef = grep {not defined $line->{$_}} ('Residue sequence number', 'Chain identifier');
		if (scalar @undef > 0) {
			p @undef;
			die "the above values are undefined from $pdb_file_with_chains line $.:\n$_";
		}
		$res_num2chain{$line->{'Residue sequence number'}} = $line->{'Chain identifier'};
	}
	close $in;
	open $in,     '<', $pdb_file_without_chains;
	open my $out, '>', $output_file;
	while (<$in>) {
		unless (/^(?:ATOM|HETATM|TER)\h/) {
			print $out $_;
			next
		}
		my $line = pdb_coord($_);
		unless (defined $res_num2chain{$line->{'Residue sequence number'}}) {
			say $_;
			p $line;
			die "line $. of $pdb_file_without_chains has an undefined atom serial number";
		}
		my @line = split '', $_, 23;
		$line[21] = $res_num2chain{$line->{'Residue sequence number'}};
		say $out join '', @line;
	}
	return 1;
}
#sub split_groups ($idx_file, $group1, $group2) {
#	my $current_sub = (split(/::/,(caller(0))[3]))[-1];
#	my $group2index = get_group_indices( $idx_file );
#	my @undef_groups = grep { not defined $group2index->{$_}} ($group1, $group2);
#	if (scalar @undef_groups > 0) {
#		p @undef_groups;
#		die "The above groups aren't defined in $current_sub";
#	}
#	foreach my $group (grep {
#										$_ ne $group1
#										&&
#										$_ ne $group2
#									} keys %{ $group2index }) {
#		delete $group2index->{$group};
#	}
#	foreach my $group ($group1, $group2) {
#		my $g_tpr = "3md.$group.tpr";
#		task({
#			cmd           => "echo $val | gmx convert-tpr -s 3md.tpr -o $g_tpr -n cpx.ndx",
#			'input.files' => ['3md.tpr', 'cpx.ndx'],
#			'log.fh'      => $log,
#			'output.files'=> $g_tpr, # only do this once
#			overwrite     => 'true'
#		});
#		my $subset_xtc = "3md.$group.$n.xtc";
#		task({
#			cmd            => "echo $val | gmx trjconv -s $g_tpr -f 3md_out$n.xtc -o $subset_xtc -n cpx.ndx",
#			'input.files'  => ["3md_out$n.xtc", $g_tpr],
#			'log.fh'       => $log,
#			'output.files' => $subset_xtc,
#			overwrite      => 'true'
#		});
#	}
#	foreach my $backup (list_regex_files('^#(?:3|chi)')) {
#		unlink $backup; # these can interfere with future commands
#		say2("deleted $backup", $log);
#	}
#}
1;
