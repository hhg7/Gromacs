#!/usr/bin/env perl

use 5.042;
no source::encoding;
use warnings FATAL => 'all';
use Test::More;
use Test::Exception; # die_ok
use Gromacs 'gromacs_log2hash';
use JSON qw(encode_json decode_json);
use List::Compare;
use Term::ANSIColor;
use Util 'dir';
use Time::HiRes;

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

sub test_hash_eq ( $logfile ) {
	my $h = gromacs_log2hash( $logfile ); 
	my $check = json_file_to_ref("$logfile.json");
	my @lkeys = keys %{ $h };
	my @rkeys = keys %{ $check };
	if (scalar @lkeys != scalar @rkeys) {
		p @lkeys;
		p @rkeys;
		die "the number of keys is not identical for $logfile";
	}
	my $lc = List::Compare->new( \@lkeys, \@rkeys);
	unless ($lc -> is_LequivalentR) {
		p @lkeys;
		p @rkeys;
		die 'the 2 sets of keys are not identical';
	}
	my @diff_keys = grep { # are all scalar valued keys equal?
							((ref $h->{$_}) eq '')
							&&
							($h->{$_} ne $check->{$_})
							} @lkeys;
	if (scalar @diff_keys > 0) {
		p @diff_keys;
		die 'the above keys are different between hash and check';
	}
	return 'equal';
}
sub list_regex_files ($regex, $directory = '.', $case_sensitive = 'yes') {
	die "\"$directory\" doesn't exist" unless -d $directory;
	my @files;
	opendir (my $dh, $directory);
	if ($case_sensitive eq 'yes') {
		$regex = qr/$regex/;
	} else {
		$regex = qr/$regex/i;
	}
	while (my $file = readdir $dh) {
		next if $file !~ $regex;
		next if $file =~ m/^\.{1,2}$/;
		my $f = "$directory/$file";
		next unless -f $f;
		if ($directory eq '.') {
			push @files, $file
		} else {
			push @files, $f
		}
	}
	@files
}
my %veritas = (
	'0em_out.log' => {
		'Finished'           => 'Finished mdrun on rank 0 Fri Dec  5 07:45:14 2025',
		'GROMACS version'		=> '2025.3',
		'Precision'				=> 'mixed',
		'Memory model'			=> '64 bit',
		'MPI library'			=> 'thread_mpi',
		'OpenMP support'		=> 'enabled (GMX_OPENMP_MAX_THREADS = 128)',
		'GPU support'			=> 'disabled',
		'SIMD instructions'	=> 'AVX2_256',
		'CPU FFT library'		=> 'fftw-3.3.8-sse2-avx',
		'GPU FFT library'		=> 'none',
		'Multi-GPU FFT'		=> 'none',
		'RDTSCP usage'			=> 'enabled',
		'TNG support'			=> 'enabled',
		'Hwloc support'		=> 'disabled',
		'Tracing support'		=> 'disabled',
		'C compiler'			=> '/usr/bin/cc GNU 11.5.0',
		'C compiler flags'	=> '-fexcess-precision=fast -funroll-all-loops -mavx2 -mfma -Wno-missing-field-initializers -O3 -DNDEBUG',
		'C++ compiler'			=> '/usr/bin/c++ GNU 11.5.0',
		'C++ compiler flags'	=> '-fexcess-precision=fast -funroll-all-loops -mavx2 -mfma -Wno-missing-field-initializers -Wno-cast-function-type-strict SHELL:-fopenmp -O3 -DNDEBUG',
		'BLAS library'			=> 'External - detected on the system',
		'LAPACK library'		=> 'External - detected on the system',
		'Host'               => 'blip-a.hpc.uidaho.edu',
		'Input Parameters' => {
			'integrator' => 'steep',
			'tinit' => 0,
			'dt' => 0.001,
			'nsteps' => 5000,
			'init-step' => 0,
			'simulation-part' => 1,
			'mts' => 'false',
			'mass-repartition-factor' => 1,
			'comm-mode' => 'Linear',
			'nstcomm' => 100,
			'bd-fric' => 0,
			'ld-seed' => -39887010,
			'emtol' => 10,
			'emstep' => 0.01,
			'niter' => 20,
			'fcstep' => 0,
			'nstcgsteep' => 1000,
			'nbfgscorr' => 10,
			'rtpi' => 0.05,
			'nstxout' => 0,
			'nstvout' => 0,
			'nstfout' => 0,
			'nstlog' => 50,
			'nstcalcenergy' => 100,
			'nstenergy' => 50,
			'nstxout-compressed' => 50,
			'compressed-x-precision' => 1000,
			'cutoff-scheme' => 'Verlet',
			'nstlist' => 10,
			'pbc' => 'xyz',
			'periodic-molecules' => 'false',
			'verlet-buffer-tolerance' => 0.005,
			'verlet-buffer-pressure-tolerance' => 0.5,
			'rlist' => 1.26,
			'coulombtype' => 'PME',
			'coulomb-modifier' => 'Potential-shift',
			'rcoulomb-switch' => 0,
			'rcoulomb' => 1.2,
			'epsilon-r' => 1,
			'epsilon-rf' => 'inf',
			'vdw-type' => 'Cut-off',
			'vdw-modifier' => 'Potential-shift',
			'rvdw-switch' => 0,
			'rvdw' => 1.2,
			'DispCorr' => 'No',
			'table-extension' => 1,
			'fourierspacing' => 0.12,
			'fourier-nx' => 80,
			'fourier-ny' => 80,
			'fourier-nz' => 80,
			'pme-order' => 4,
			'ewald-rtol' => 1e-05,
			'ewald-rtol-lj' => 0.001,
			'lj-pme-comb-rule' => 'Geometric',
			'ewald-geometry' => '3d',
			'epsilon-surface' => 0,
			'ensemble-temperature-setting' => 'not available',
			'tcoupl' => 'No',
			'nsttcouple' => -1,
			'nh-chain-length' => 0,
			'print-nose-hoover-chain-variables' => 'false',
			'pcoupl' => 'No',
			'refcoord-scaling' => 'No',
		}
	},
	'3md_out03.log' => {
		Host => 'n123.fortytwo.ibest.uidaho.edu',
		'Input Parameters' => {
			'integrator' => 'md',
			'tinit' => 0,
			'dt' => 0.002,
			'nsteps' => 500000000,
			'init-step' => 0,
			'simulation-part' => 1,
			'mts' => 'false',
			'mass-repartition-factor' => 1,
			'comm-mode' => 'Linear',
			'nstcomm' => 100,
			'bd-fric' => 0,
			'ld-seed' => -1511134482,
			'emtol' => 10,
			'emstep' => 0.01,
			'niter' => 20,
			'fcstep' => 0,
			'nstcgsteep' => 1000,
			'nbfgscorr' => 10,
			'rtpi' => 0.05,
			'nstxout' => 0,
			'nstvout' => 0,
			'nstfout' => 0,
			'nstlog' => 500,
			'nstcalcenergy' => 100,
			'nstenergy' => 500,
			'nstxout-compressed' => 500,
			'compressed-x-precision' => 1000,
			'cutoff-scheme' => 'Verlet',
			'nstlist' => 40,
			'pbc' => 'xyz',
			'periodic-molecules' => 'false',
			'verlet-buffer-tolerance' => 0.005,
			'verlet-buffer-pressure-tolerance' => 0.5,
			'rlist' => 1.281,
			'coulombtype' => 'PME',
			'coulomb-modifier' => 'Potential-shift',
			'rcoulomb-switch' => 0,
			'rcoulomb' => 1.2,
			'epsilon-r' => 1,
			'epsilon-rf' => 'inf',
			'vdw-type' => 'Cut-off',
			'vdw-modifier' => 'Potential-shift',
			'rvdw-switch' => 0,
			'rvdw' => 1.2,
			'DispCorr' => 'EnerPres',
			'table-extension' => 1,
			'fourierspacing' => 0.12,
			'fourier-nx' => 80,
			'fourier-ny' => 80,
			'fourier-nz' => 80,
			'pme-order' => 4,
			'ewald-rtol' => 1e-05,
			'ewald-rtol-lj' => 0.001,
			'lj-pme-comb-rule' => 'Geometric',
			'ewald-geometry' => '3d',
			'epsilon-surface' => 0,
			'ensemble-temperature-setting' => 'constant',
			'ensemble-temperature' => 300,
			'tcoupl' => 'V-rescale',
			'nsttcouple' => 100,
			'nh-chain-length' => 0,
			'print-nose-hoover-chain-variables' => 'false',
			'pcoupl' => 'C-rescale',
			'pcoupltype' => 'Isotropic',
			'nstpcouple' => 100,
			'tau-p' => 5,
		}
	}
);
foreach my $test_log ('0em_out.log', '3md_out03.log','md.log', '3md_out04.log') {
	my $t0 = Time::HiRes::time();
	my $h = gromacs_log2hash($test_log);
	my $t1 = Time::HiRes::time();
	printf("Did $test_log in %lf seconds\n", $t1-$t0);
#	p $h->{colvars} if defined $h->{colvars};
	next unless defined $veritas{$test_log};
	my @err = grep { # missing keys
							not defined $h->{$_}
						} keys %{ $veritas{$test_log} };
	if (scalar @err == 0) {
		ok(1, "$test_log isn't missing any defined truth keys");
	} else {
		p @err;
		die "The above keys are missing from the hash made from $test_log";
	}
	@err = grep { # checking values of scalar keys
						(ref $veritas{$test_log}{$_} eq '')
						&&
						($veritas{$test_log}{$_} ne $h->{$_})
						} keys %{ $veritas{$test_log} };
	if (scalar @err == 0) {
		ok(1, "$test_log scalar keys are all correct");
	} else {
		p @err;
		my @val = @{ $veritas{$test_log} }{@err};
		say STDERR "These are the \"true\" values for $test_log:";
		p @val;
		@val = @{ $h }{@err};
		say STDERR "But these values were detected from $test_log:";
		p @val;
		die "The above keys are incorrect in $test_log";
	}
	foreach my $key ('Input Parameters') {
		@err = grep { # missing keys
							not defined $h->{$key}{$_}
						} keys %{ $veritas{$test_log}{$key} };
		if (scalar @err == 0) {
			ok(1, "$test_log has all scalar keys for \"$key\"");
		} else {
			p @err;
			die "The above keys are missing from the hash made from $test_log and \"$key\"";
		}
		@err = grep { # checking values of scalar keys
						(ref $veritas{$test_log}{$key}{$_} eq '')
						&&
						($veritas{$test_log}{$key}{$_} ne $h->{$key}{$_})
						} keys %{ $veritas{$test_log}{$key} };
		if (scalar @err == 0) {
			ok(1, "$test_log -> \"$key\" scalar keys are all correct");
		} else {
			p @err;
			die "The above keys are incorrect in $test_log for \"$key\"";
		}
	}
}
done_testing();
