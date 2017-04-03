#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Cwd 'abs_path';
use Bio::TreeIO;

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
my @genome = @ARGV;
my $allarg = join(' ', @genome);
prodigal_multi($allarg);

## Run HMMs, align, assemble ML
my %famseq = ();
foreach my $faa (@genome){
    $faa =~ s/\.fna$/\.faa/;
    my $fsr = runhmm($faa, \%famseq);
    %famseq = %$fsr;
}
alignall(\%famseq);
my %multilocus = ();
foreach my $tigr (glob("TIGR*.fs.afa")){
    my $afa = new Bio::SeqIO(-file=>$tigr, -format=>'fasta');
    my %tmp = ();
    while(my $seq=$afa->next_seq){
	my ($gene, $id) = split(/--/, $seq->id);
	$id =~ s/.+\/(.+)\.faa$/$1/;
	$multilocus{$id} .= $seq->seq;
    }
}
open my $mfh, '>', 'multilocus.afa' or die $!;
foreach my $id (keys %multilocus){
    print $mfh '>'.$id."\n".$multilocus{$id}."\n";
}
close $mfh;

## Tree
print "Generating tree...";
system("FastTree -log multilocus.log -quiet < multilocus.afa > multilocus.nwk");
print "DONE!\n";

## Cleanup
#my @torm = ('tmp*', '*.reduced', 'RAxML_*', 'bs-files', '*.fs.*', '*.prod.*', 'allfam.tre');
#foreach (@torm){
    #system("rm $_");
#}

sub prodigal_multi{
    my $big_list = shift;
    print STDERR "Calling genes with prodigal...";
    system("python $script_dir/prodigal_multi.py $big_list");
    print STDERR "DONE!\n";
}

sub runhmm{
    my ($faa, $fsr) = (shift, shift);
    my %fs = %$fsr;
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
    foreach my $tf (sort keys %cut){
	if(exists $ann{$tf}{'q'}){
	    $fs{$tf}{$pref}{'seq'} = $s{$ann{$tf}{'q'}};
	    $fs{$tf}{$pref}{'orf'} = $ann{$tf}{'q'};
	}else{
	    $fs{$tf}{$pref}{'orf'} = 'none';
	    foreach(1..$hmmlen{$tf}){
		$fs{$tf}{$pref}{'seq'} .= 'X';
	    }
	}
    }
    print "DONE!\n";
    return(\%fs);
}

sub alignall{
    my $fsr = shift;
    my %fs = %$fsr;
    my $f = scalar(keys %fs);
    foreach my $t (keys %fs){
	print STDERR "$t\tAligning individual gene family...";
	open my $fsh, '>', "$t.fs.faa" or die $!;
	foreach my $g (sort keys %{$fs{$t}}){
	    print $fsh '>' . join('--', $fs{$t}{$g}{'orf'}, $g) . "\n" . $fs{$t}{$g}{'seq'} . "\n";
	}
	close $fsh;
	system("mafft --quiet --namelength 70 $t.fs.faa > $t.fs.afa");
	print "DONE!\n";
    }
}
