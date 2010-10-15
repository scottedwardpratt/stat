#!/bin/sh -- # A comment to find the right perl
eval 'exec perl -S $0 ${1+$@}'
    if 0;

use Getopt::Long;

sub func {
  # parameter definitions
  # p[0] = 1/T^2
  # p[1] = 1/T^4
  # p[2] = Tanh offset
  # p[3] = Tanh slope
  # p[4] = Tanh power
    my $x = shift @_;
    my @p = @_;
    my $f;

    # print "@p\n"; 
    my $h = ($p[0]/($x**2))+($p[1]/($x**4));  # hi-T part
    my $g = 1-1/(1+exp(($x-$p[2])/$p[3]))**$p[4];    # tanh suppression
    my $f = $g*$h;
    
   return $f;
}

GetOptions ( "help"      => \$help
	     );

if (($help)||(@ARGV !=4)) {
    print(
	  "Script by M. Cheng to convert eos parameterizations into output suitable for VH2 hydro.\n",
	  "Writes out Temperature, Energy, and Pressure (/ T^4).\n",
	  "Usage: $0 [options] paramfile T_0(GeV) step_size(GeV) T_Max (GeV)\n",
	  "Options:\n",
	  "        --help       Print this message.\n",
	  );
    exit;
}

my $paramfile = shift @ARGV;
my $Tstart = shift @ARGV;
my $interval = shift @ARGV;
my $Tmax = shift @ARGV;

open FILE, "<$paramfile" or die "Could not open $paramfile for input:$!\n";
while (<FILE>) {
    chomp;
    my @l = split;
    # print "@l\n";
    my $tag = shift @l;
    @{$params{$tag}} = @l;
}
close FILE;

foreach $tag (keys %params) {
  $p{$tag}[0] = 0;
  $int{$tag}[0] = &func($Tstart,@{$params{$tag}});
  $e{$tag}[0] = $int{$tag}[0];
  $T = $Tstart;
  $tfirst = 1;

  $numtemps = int(($Tmax - $T)/$interval) + $tfirst;
  for ($t = $tfirst; $t < $numtemps+2; $t++) {
    $T += $interval;
    $T{$tag}[$t] = $T;
    $int{$tag}[$t] = &func($T,@{$params{$tag}});
    my $aux = &func($T-$interval/2,@{$params{$tag}});
    $p{$tag}[$t] = $p{$tag}[$t-1] + ($interval/6)*($int{$tag}[$t-1]+4*$aux+$int{$tag}[$t])/$T;
    $e{$tag}[$t] = $int{$tag}[$t]+3*$p{$tag}[$t];
    $s{$tag}[$t] = $e{$tag}[$t]+$p{$tag}[$t];
    $p_over_e{$tag}[$t] = $p{$tag}[$t]/$e{$tag}[$t];

#   Calculate cs2 for $t-1 (this will overwrite previous entry for merger)
    $cs2{$tag}[$t-1] = ($p{$tag}[$t]*($T)**4-$p{$tag}[$t-2]*($T-(2*$interval))**4)/($e{$tag}[$t]*($T)**4-$e{$tag}[$t-2]*($T-(2*$interval))**4);

  }
}

for $tag (keys %params) {
    my $outputfile = "$tag"."_eos.dat";
#    print "Tag = $tag\n";
    open OUT, ">$outputfile" or die "Could not open $outputfile for output:$!\n";
    print OUT "#T Interaction energy_density pressure entropy p_over_e cs^2\n";
#    print OUT "$numtemps\n";
    for ($t = 1; $t < $numtemps+1; $t++) {
	print OUT "$T{$tag}[$t] $int{$tag}[$t] $e{$tag}[$t] $p{$tag}[$t] $s{$tag}[$t] $p_over_e{$tag}[$t] $cs2{$tag}[$t] \n";
#	print "$tag : $T{$tag}[$t] $e{$tag}[$t] $p{$tag}[$t] $cs2{$tag}[$t]\n";
#	print OUT "$T{$tag}[$t] $e{$tag}[$t] $p{$tag}[$t] $cs2{$tag}[$t]\n";
    }
    close OUT;

}
