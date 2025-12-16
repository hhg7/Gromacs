require 5.036;
use 5.042;
no source::encoding;
use DDP { output => 'STDOUT', array_max => 10, show_memsize => 1 };
use Devel::Confess 'color';
package Gromacs; # packages after this are used, but aren't available when the module is loaded
our $VERSION = 0.01;
use File::Temp 'tempfile';
use List::MoreUtils 'first_index';
use DDP { output => 'STDOUT', array_max => 10, show_memsize => 1 };
use Devel::Confess 'color';
use Cwd 'getcwd';
use Exporter qw(import);
#use Term::ANSIColor;
our @EXPORT = qw(gromacs_log2hash plot_gromacs_energy make_rmsd_json get_group_indices);
our @EXPORT_OK = @EXPORT;
use JSON 'encode_json';
use Matplotlib::Simple;
use SimpleFlow;

sub gromacs_log2hash ($log, $json_filename = undef) {
	open my $fh, '<', $log;
	my @log = <$fh>;
	close $fh;
	splice @log, 0, first_index { /^GROMACS:\h+/ } @log;
	my %data;
	@log = grep {$_ !~ m/^\++\h+PLEASE .+\++/} @log;
	@log = grep {$_ !~ m/^DOI:\h/} @log;
	@log = grep {$_ !~ m/^\-+.+Thank You\h\-+/} @log;
	foreach my $i (grep {$log[$_] =~ m/^DD\h+step\h+\d+/} reverse 0..$#log) {
		splice @log, $i-1, 2; # remove this line and the blank that comes before it
	}
	foreach my $i (grep {$log[$_] =~ m/^step\h+\d+\h+Turning on/} reverse 0..$#log) {
		splice @log, $i, 1;
	}
	foreach my $i (grep {$log[$_] =~ m/^Current ref_t for group System:\h+\d+\.\d+/} reverse 0..$#log) {
		splice @log, $i-1, 2;
	}
	@log = grep {$_ !~ m/^colvars:\h+Initializing\h/} @log;
	foreach my $i (reverse 0..$#log) {
		if ($log[$i] =~ m/^\h+Time:\h+(\d+)\.(\d+)\h+(\d+)\.(\d+)\h+(\d+)\.(\d+)/) {
			$data{'Core Time (s)'} = "$1.$2";
			$data{'Wall Time (s)'} = "$3.$4";
			$data{'Time (%)'} = "$5.$6";
			splice @log, $i, 1;
			$data{'Wall Time (min)'} = $data{'Wall Time (s)'} / 60;
			$data{'Wall Time (hr)'}  = $data{'Wall Time (s)'} / 3600;
			last;
		}
	}
	my $str = join ('я', @log);
	my $performance_i = first_index { /^Performance:\h+\d+/ } @log;
	if (
			($performance_i > 0)
			&&
			($log[$performance_i] =~ m/^Performance:\h+(\d+)\.(\d+)\h+(\d+)\.(\d+)/)
		) {
		$data{'ns/day'}  = "$1.$2";
		$data{'hour/ns'} = "$3.$4";
	}
	my $last_i = first_index { /^\h+Statistics over \d+ steps using \d+ frames/ } @log;
	if ($last_i > 0) {
		splice @log, $last_i - $#log;
	}
	foreach my $line (grep {/^There are: \d+ Atoms$/} @log) {
		if ($line =~ m/(\d+)/) {
			$data{'Atom Count'} = $1;
		} else {
			die "$line failed regex.";
		}
	}
	chomp @log;
	foreach my $line (grep {$_ =~ m/^.+:\h+\H+/} @log) {
		my @line = split /:\h+/, $line;
		$data{$line[0]} = $line[1];
	}
#	p @log, array_max => scalar @log;
	foreach my $writing_index (reverse grep { $log[$_] =~ m/^Writing checkpoint/} 0..$#log) {
		splice @log, $writing_index, 3;
	}
	my $input_param_i = first_index {$_ eq 'Input Parameters:'} @log;
	my $qm_opts_i     = first_index {$_ eq 'qm-opts:'}          @log;
	foreach my $i (reverse $input_param_i..$qm_opts_i) {
		$log[$i] =~ s/^\h+//;
		my @line = split /\h+=\h+/, $log[$i];
		splice @log, $i, 1; # prevent regex confusion later on
		next unless scalar @line == 2;
		$data{$line[0]} = $line[1];
	}
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
									
								} 0..$#log;
	if (scalar @time_indices == 0) {
		p @log, array_max => scalar @log;
		my @step_time_i = grep {$log[$_] =~ m/^\h+Step\h+Time$/} 0..$#log;
		my @t = @log[@step_time_i];
		p @t;
		die "Couldn't get times for $log";
	}
	my @energies;
	foreach my $time_index (@time_indices) {
		if ($log[$time_index] =~ m/(\d+)\.(\d+)$/) {
			push @{ $data{time} }, "$1.$2";
		} else {
			die "$log[$time_index] failed regex.";
		}
	}
	my $reading_energies = 'false';
	foreach (@log) {
		if ($_ eq '   Energies (kJ/mol)') {
			$reading_energies = 'true';
			next;
		}
		if ($_ eq '') {
			$reading_energies = 'false';
			next;
		}
		next unless $reading_energies eq 'true';
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
	foreach my $i (grep {$log[$_] =~ m/^Running on \d+ node/} reverse 0..$#log) {
		if ($log[$i] =~ m/^Running on (\d+) nodes? with total (\d+) cores, (\d+) processing units$/) {
			$data{nodes} = $1;
			$data{cores} = $2;
			$data{'processing units'} = $3;
		}
		splice @log, $i, 1;
	}
	foreach my $i (grep {$log[$_] =~ m/^Hardware detected on host /} reverse 0..$#log) {
		if ($log[$i] =~ m/^Hardware detected on host (.+):$/) {
			$data{Host} = $1;
		}
		splice @log, $i, 1; # don't possibly confuse future regex
	}
	if ($str =~ m/
	\h+Vendor:\h+[^я]+я
	\h+Brand:\h+([^я]+)
	/x) {
		$data{CPU} = $1;
	} else {
		p $str;
		die "$log failed CPU info regex.";
	}
	foreach my $key (grep {/^\h+/} keys %data) {
		my $new_key = $key;
		$new_key =~ s/^\h+//;
		$data{$new_key} = delete $data{$key};
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
#and returns
#    [0] "System",
#    [1] "Protein",
#    [2] "Receptor",
#    [3] "Ligand_CA",
#    [4] "Ligand",
#    [5] "BindingSite_CA",
#    [6] "BindingSite",
#    [7] "Receptor-BindingSite",
#    [8] "BindingAxisAnchor"
	my $t = task({
		cmd           => "gmx check -n $ndx_file",
		'input.files' => $ndx_file,
		note          => 'read ndx file'
	});
	my @ndx = split /\n/, $t->{stdout};
	my $i = first_index {/^\h+\d+\h/} @ndx;
	splice @ndx, 0, $i;
	my @i2group;
	foreach my $line (@ndx) {
		my @line = grep {$_ =~ m/\H/} split /\h+/, $line, 4;
		push @i2group, $line[1];
	}
	p @i2group, show_memsize => 1;
	return \@i2group;
}
1;
