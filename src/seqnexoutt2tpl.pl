#!/usr/local/bin/perl
use strict;


# option process (help message)

if ($ARGV[0] eq "-h") {
	print <<EOF;
This program convert xxx.seq.nex.out.t to xxx.tpl
    Usage: cat xxx.seq.nex.out.t | seqnexoutt2tpl.pl [option] > xxx.tpl
    Option: -h    show this help message
            -c    cut brunch length information
EOF

    exit;
}

my $option = 0;
if ($ARGV[0] eq "-c") {
	$option = 1;
}


# header process

my @translate_index = ();

while (<STDIN>) {
	chomp();
	if (/^[\s]*translate[\s]*$/) {
		last;
	}
}

while (<STDIN>) {
	chomp();
	if (/(\d+)\s+(\w+)/) {
		#print "$1, $2\n";
		$translate_index[$1] = $2;
	}
	if (/;/) {
		last;
	}
}


# data input from stdin and filtering

my @tree_data = ();
my $tree_count = 0;

while (<STDIN>) {
#for (1..1) {
#	$_ = <STDIN>;
	chomp();
	if (/^end;$/) {
		last;
	}
	s/^[^\(]+//g;
    #print $_, "\n";
	s/(\d+):/$translate_index[$1]:/g;
	my $tree = $_;
	if ($option == 1) {
		$tree =~ s/:[\d\.]+//g;
	}
	#print $tree, "\n";
	push(@tree_data, $tree);
	$tree_count ++;

}


# data output to stdout

print $tree_count, "\n";
$tree_count = 1;
my $i;
foreach $i (@tree_data) {
	print $i, " ", $tree_count, "\n";
	$tree_count ++;
}

exit;
