#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Cwd 'abs_path';
use Bio::TreeIO;

my ($pref, $genome, $refpkg, $refaln) = (shift,shift,shift,shift);

## Set params
our $cpu = 8;

## Setup script path
my $script_dir = abs_path($0);
my @p = split(/\//, $script_dir);
$script_dir =~ s/$p[-1]//;

## Set paths
our $tigrdir = "$script_dir/tigrfam";
our $tigrdb = $tigrdir . '/genprop0799.hmmdb';
our $tigrcut = $tigrdir . '/genprop0799.cutoffs.tsv';

## Read in all cutoffs
my %cut = ();
open my $tch, '<', $tigrcut or die $!;
while(<$tch>){
	unless($_ =~ m/^#/){
		my ($model, $tc, $nc) = split(/\t/, $_);
		$cut{$model} = {
			'trust'	=> $tc,
			'noise'	=> $nc
		};
	}
}
close $tch;

## Read in all hmm lengths
my %hmmlen = ();
open my $tdbh, '<', $tigrdb or die $!;
my $cur = '';
while(<$tdbh>){
	chomp;
	if($_ =~ m/^NAME\s+(\S+)/){
		$cur = $1;
	}elsif($_ =~ m/^LENG\s+(\d+)/){
		$hmmlen{$cur} = $1;
	}
}
close $tdbh;

## Call genes
prodigal_multi($genome);

## Run HMMs, align, assemble ML
my $faa = $genome;
$faa =~ s/\.fna$/\.faa/;
my $multilocus = runhmm($faa);
open my $mfh, '>', "$pref.ml.faa" or die $!;
print $mfh '>'.$pref."\n".$multilocus."\n";
close $mfh;

## Align to ref
print "Generating alignment...";
system("cat $pref.ml.faa $refaln > $pref.ref.faa");
system("mafft --quiet --namelength 70 $pref.ref.faa > $pref.ref.aln.fasta");
print "DONE!\n";

## Place and get tree
print "Placing on tree...";
system("pplacer -c $refpkg $pref.ref.aln.fasta");
system("guppy sing $pref.ref.aln.jplace");
print "DONE!\n";

sub prodigal_multi{
    my $big_list = shift;
    print STDERR "Calling genes with prodigal...";
    system("python $script_dir/prodigal_multi.py $big_list");
    print STDERR "DONE!\n";
}

sub runhmm{
    my $faa = shift;
    print STDERR "$faa\tScanning with TIGRFAM HMMs...";
    ## HMMscan of all prodigal faas
    my $pref = $faa;
    $pref =~ s/\.prod\.faa//;
    system("hmmscan -o tmp.hmmscan.out --tblout tmp.hmmtbl.out -E 1e-5 --cpu $cpu --noali $tigrdb $faa"); ## COMMENT OUT TO SKIP
    ## Parse hmmtbl
    open my $sfh, '<', "tmp.hmmtbl.out" or die "Died in hmmscanner: $!";
    my %ann = ();
    my %bit = ();
    while(<$sfh>){
	unless($_ =~ m/^#/ || $_ =~ m/^\W/){
	    chomp;
	    my ($hit, $mod, $query, $dash, $evalue, $score, @rest) = split(/\s+/, $_);
	    $bit{$query}{$hit} = $score;
	}
    }
    close $sfh;
    ## Assign annotations
    my %grab = ();
    foreach my $q (sort keys %bit){
	foreach my $h (sort keys %{$bit{$q}}){
	    if( $bit{$q}{$h} >= $cut{$h}{'trust'}){
		if(exists $ann{$h}{'q'}){
		    if($ann{$h}{'b'} < $bit{$q}{$h}){
			$ann{$h}{'q'} = $q;
			$ann{$h}{'b'} = $bit{$q}{$h};
			$grab{$q} = 1;
		    }
		}else{
		    $ann{$h}{'q'} = $q;
		    $ann{$h}{'b'} = $bit{$q}{$h};
		    $grab{$q} = 1;
		}
	    }
	}
    }
    ## Grab protein seqs
    my %s = ();
    my $fa = new Bio::SeqIO(-file=>$faa, -format=>'fasta');
    while(my $seq = $fa->next_seq){
	if(exists $grab{$seq->id}){
	    my $ss = $seq->seq;
	    $ss =~ s/\*//;
	    $s{$seq->id} = $ss;
	}
    }
    my $ml = '';
    foreach my $tf (sort keys %cut){
	if(exists $ann{$tf}{'q'}){
	    $ml .= $s{$ann{$tf}{'q'}};
	}else{
	    foreach(1..$hmmlen{$tf}){
		$ml .= 'X';
	    }
	}
    }
    print "DONE!\n";
    return($ml);
}
