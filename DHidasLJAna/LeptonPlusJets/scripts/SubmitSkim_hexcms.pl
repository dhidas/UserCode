#!/usr/bin/perl -w

use File::Copy;

if ($#ARGV != 1) {
  print "Usage: [INDIR] [OUTDIR]\n";
  exit;
}

my $INDIR = $ARGV[0];
my $OUTDIR = $ARGV[1];


my $NAME = "Skim";
my $LISTNAME = "$OUTDIR/$NAME.list";

my $RELEASEDIR = `pwd`;
$RELEASEDIR =~ m/(.*)\/src\//;
$RELEASEDIR = "$1/src";


print "NAME:       $NAME\n";
print "INDIR:      $INDIR\n";
print "OUTDIR:     $OUTDIR\n";
print "RELEASEDIR: $RELEASEDIR\n";


mkdir "$OUTDIR" or die "$OUTDIR already exists.  Move it first and then try again.\n";
mkdir "$OUTDIR/Log" or die "cannot make out/Log dir $!";

my $EXENAME = "$RELEASEDIR/DHidasLJAna/LeptonPlusJets/scripts/run_hexcms.sh";
my $PYFILE  = "$RELEASEDIR/DHidasLJAna/LeptonPlusJets/python/run/Skim_Template_cfg.py";
my $SRCFILE = "$RELEASEDIR/DHidasLJAna/LeptonPlusJets/plugins/DHidasPatAna.cc";
copy $PYFILE,  $OUTDIR;
copy $SRCFILE, $OUTDIR;

my @FILES = `ls -1 $INDIR/*.root`;
#my @FILES = `cat ZPrime.list`;
foreach (@FILES) {
  print "Using file: $_";
}

my $NFiles = $#FILES + 1;
print "Number of files found: $NFiles\n";
my $NSECTIONS = 12;
print "Number of sections to queue: $NSECTIONS\n";


if (-e $LISTNAME) {
  print "$LISTNAME exists.  Do you want to overwrite [y/n]: ";
  my $YN = <STDIN>;
  if ($YN =~ m/^y/i) {
  } else {
    exit;
  }
}

open LIST, ">$LISTNAME" or die "cannot open the list file $LISTNAME $!";
print LIST @FILES;
close LIST;

my $SUBMITNAME = "$OUTDIR/Submit_$NAME.jdl";
open SUBMIT, ">$SUBMITNAME" or die "cannot open the submit file for writing $!";
print SUBMIT "universe = vanilla\n";
print SUBMIT "+AccountingGroup = \"group_rutgers.dhidas\"\n";
print SUBMIT "Executable = $EXENAME\n";
print SUBMIT "\n";
print SUBMIT "Should_Transfer_Files = NO\n";
print SUBMIT "Output = $OUTDIR/Log/Log_\$(Process).stdout\n";
print SUBMIT "Error =  $OUTDIR/Log/Log_\$(Process).stderr\n";
print SUBMIT "Log =    $OUTDIR/Log/Log_\$(Process).log\n";
print SUBMIT "notify_user = dhidas\@physics.rutgers.edu\n";
print SUBMIT "Arguments = \$(Process) $OUTDIR $LISTNAME $NESCTIONS $RELEASEDIR\n";
print SUBMIT "Queue $NSECTIONS\n";
close SUBMIT;


print "condor_submit $SUBMITNAME\n";
`condor_submit $SUBMITNAME`;
