#!/usr/bin/env perl

use 5.042.1;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Util;
#use latex qw(write_2d_array_to_tex_tabular);
# s/[!@#\$\%^&*\(\)\{\}\[\]\<\>,\/'"\-\h;\+=]+/_/g; # annoying chars
my (@undef, %diff);
foreach my $logfile (dir(regex => 'md.+log$', dir => "$ENV{HOME}/ui/pepPriML/", recursive => 1)) {
	my $total_lines = 0;
	my $finished_line;
	open my $fh, '<', $logfile;
	while (<$fh>) {
		$total_lines++;
		if (/^Finished\h/) {
			$finished_line = $.;
		}
	}
	if (defined $finished_line) {
		$diff{ $total_lines - $finished_line } ++;
	} else {
		push @undef, $logfile;
	}
}
p %diff;
p @undef;
