#!/usr/bin/env perl

use 5.042;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Devel::Confess 'color';
use Scalar::Util 'looks_like_number';

open my $fh, '<', '3md_out03.log';
while (<$fh>) {
	last if $. > 240;
	next unless (174 <= $. <= 239);
	chomp;
	$_ =~ s/^\h+//;
	my @line = split /\h+=\h+/, $_;
	next unless scalar @line == 2;
	if (looks_like_number($line[1])) {
		say "'$line[0]' => $line[1],";
	} else {
		say "'$line[0]' => '$line[1]',";
	}
}
