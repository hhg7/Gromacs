#!/usr/bin/env perl

use 5.042;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':all';
use Devel::Confess 'color';
use Time::HiRes;
use List::MoreUtils 'lastidx';

open my $fh, '<', '3md_out03.log';
my @log = <$fh>;
close $fh;
my $t0 = Time::HiRes::time();
for (my $n = 0; $n < 999; $n++) {
	my $i = lastidx { /^Performance:\h+\d+/ } @log;
}
my $t1 = Time::HiRes::time();
printf("lastidx = %g\n", $t1-$t0);
$t0 = Time::HiRes::time();
for (my $n = 0; $n < 999; $n++) {
	foreach my $i (reverse 0..$#log) {
		last if $log[$i] =~ m/^Performance:\h+\d+/;
	}
}
$t1 = Time::HiRes::time();
printf("reverse on foreach/grep = %g\n", $t1-$t0);
