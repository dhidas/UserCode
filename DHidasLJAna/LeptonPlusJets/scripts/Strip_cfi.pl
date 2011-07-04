#!/usr/bin/perl -w

# This file is used to look at a cfi file with a list of input files
# and strip it to just a text list of files

# Argument is filename you wish to strip
my $INFILE = $ARGV[0];
#print "INFILE: $INFILE\n";

# Open the input file
open IN, "<$INFILE" || die "ERROR: cannot open input file: $INFILE";

# Loop over all lines looking for a match
foreach (<IN>) {
  if (m/file:(.*)\.root/) {
    print "$1.root\n";
  }
}

