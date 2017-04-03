#!/bin/env perl

use strict;
use warnings;

my ($pref, $tree, $aln, $log) = (shift, shift, shift, shift);
system("taxit create -l $pref -P $pref.refpkg --aln-fasta $aln --tree-stats $log --tree-file $tree");
