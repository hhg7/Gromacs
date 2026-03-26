#!/usr/bin/env perl

use 5.042.1;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Gromacs 'gromacs_log2hash';

if (not defined $ARGV[0]) {
	die 'no argument was provided to ' . __FILE__;
}
my $h = gromacs_log2hash( $ARGV[0] );
if (defined $ARGV[1]) {
	p $h->{$ARGV[1]};
} else {
	p $h;
}
