#!/usr/bin/perl -w

if ($#ARGV ne 1) {
  print "Usage: ./MakeFilePY.pl [InDir] [OutFileName]\n";
  exit;
}

open OUTFILE, ">$ARGV[1]" || die;

@FI = `/bin/ls -1 $ARGV[0]/*.root`;
$Last = $FI[$#FI];
chomp $Last;

print OUTFILE "import FWCore.ParameterSet.Config as cms\n";
print OUTFILE "\n\n\n";
print OUTFILE "source = cms.Source(\"PoolSource\",\n";
print OUTFILE "                            fileNames = cms.untracked.vstring(\n";


foreach (@FI) {
  chomp;
  if ($_ ne $Last) {
    print OUTFILE "                             'file:$_',\n"
  } else {
    print OUTFILE "                             'file:$_'\n"
  }
}

print OUTFILE "                                                             )";
print OUTFILE "                              )";


print "Last File index: $#FI\n";

close OUTFILE;

