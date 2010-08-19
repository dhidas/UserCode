#!/usr/bin/perl


my $DIR = $ARGV[0];
print "Directory: $DIR\n";

my @EPS = `ls $DIR/\*.eps`;

foreach (@EPS) {
  chomp;
  print "\\includegraphics[width=0.5\\linewidth]{$_}\n";
}



