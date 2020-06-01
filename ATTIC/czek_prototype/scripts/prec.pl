#!/usr/bin/perl

$line = <>;
while($line = <>) {
$line =~ s/\r//g;
@inf = split /[\t]/, $line;
$id = shift @inf; print "$id";
foreach(@inf) { printf("\t%.2f", $_); }
print "\n";
}
