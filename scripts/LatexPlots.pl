#!/usr/bin/perl -w

my $TEMPFILENAME = "Plots";

my @EPS;



if (-e $TEMPFILENAME.".tex") {
  die "$TEMPFILENAME already exists.  move it first";
}

open MYFILE, ">$TEMPFILENAME.tex";

print MYFILE "\\documentclass[prl]{revtex4}\n";
print MYFILE "\\usepackage{graphicx}\n";
print MYFILE "\\begin{document}\n\n";

foreach (@ARGV) {
  print "Directory: $_\n";
  @EPS = `ls $_/\*.eps`;

  s/_/\\_/g;
  print MYFILE "\\section{$_}\n";
  foreach (@EPS) {
    chomp;
    print MYFILE "\\includegraphics[width=0.5\\linewidth]{$_}\n";
  }
}



print MYFILE "\\end{document}";
close MYFILE;

`latex $TEMPFILENAME.tex`;
`latex $TEMPFILENAME.tex`;
`dvips -t letter -Ppdf -G0 -tletter $TEMPFILENAME.dvi -o $TEMPFILENAME.ps`;
`ps2pdf $TEMPFILENAME.ps`;
`rm $TEMPFILENAME.ps $TEMPFILENAME.dvi $TEMPFILENAME.tex $TEMPFILENAME.log $TEMPFILENAME.aux`;


