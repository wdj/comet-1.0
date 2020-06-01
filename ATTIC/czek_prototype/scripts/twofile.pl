#!/usr/bin/perl

open FH, "foo"or die("bleh\n");
open GH, "test.good"or die("bleh\n");
$ctr = 0;
while($line = <FH>) {
  $line2 = <GH>;
  @inf1 = split /[\t\n]+/, $line;
  @inf2 = split /[\t\n]+/, $line2;
  $ctr++;
  for($i = 0; $i < @inf1; $i++) {
    if($inf1[$i] ne $inf2[$i]) { print "Line $ctr, field $i: $inf1[$i] $inf2[$i]\n"; }
  }
}
close FH;
close GH;
