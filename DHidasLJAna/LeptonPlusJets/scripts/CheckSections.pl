#!/usr/bin/perl -w


my @stuff = ();
foreach (@ARGV) {
  s/...*_(\d+)_(\d+)_...\.root/$1/;
  push(@stuff,$_);
}

@stuff = sort @stuff;

my $last = -1;
foreach (@stuff) {
  if ($_ eq $last) {
    print "Duplicate section: $_\n";
  }
  $last = $_;
}


