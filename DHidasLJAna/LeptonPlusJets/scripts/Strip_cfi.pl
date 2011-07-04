#!/usr/bin/perl -w

my $INFILE = $ARGV[0];

print "INFILE: $INFILE\n";


open IN, "<$INFILE" || die "ERROR: cannot open input file: $INFILE";

foreach (<IN>) {
  if (m/file:(.*)\.root/) {
    print "$1.root\n";
  }
}

