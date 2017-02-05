#!/usr/bin/perl
#04 February 2017 - Adam D Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

my $usage = 'perl crossReferenceClinVar.pl <genome> <clinvar> <output> 
';

die $usage , unless @ARGV == 3;
my ( $genome , $clinvar , $output ) = @ARGV;

my $IN1 = FileHandle->new( "$genome" , "r" );
if ( not defined $IN1 ) { die "ADSERROR: Could not open/read $genome\n"; }

my $IN2 = FileHandle->new( "$clinvar" , "r" );
if ( not defined $IN2 ) { die "ADSERROR: Could not open/read $clinvar\n"; }

my $OUT = FileHandle->new( "$output" , "w" );
if ( not defined $OUT ) { die "ADSERROR: Could not open/write $output\n"; }

my %snps;
my %rsIDs;
while ( my $line = $IN1->getline ) {
	next if ( $line =~ m/^#/ );
	chomp( $line );
	my ( $rsID , $chr , $pos , $genotype ) = split( "\t" , $line );
	my ( $allele1 , $allele2 ) = split // , $genotype;
	if ( not defined $allele2 ) { $allele2 = "."; }
	my $gen = join( ":" , ( $chr , $pos ) );
 	$snps{$gen} = $allele1.":".$allele2;
	$rsIDs{$gen} = $rsID;
}
$IN1->close();

$OUT->print( "RefStatus\tAltStatus\trsID\tChromosome\tPosition\tAllele1\tAllele2\t".$IN2->getline."\n" );
while ( my $line = $IN2->getline ) {
	chomp( $line );
	my ( $chr , $pos , $ref , $alt ) = split( "\t" , $line );
	if ( $line =~ /Pathogenic/ig and $line =~ /Benign/ig ) {
		$status = "C";
	}
	my $gen = join( ":" , ( $chr , $pos ) );
	if ( exists $snps{$gen} ) {
		my ( $refCount , $altCount , $rstatus , $astatus ) = (0)x4;
		my ( $allele1 , $allele2 ) = split( ":" , $snps{$gen} );
		$rstatus = "missing";
		if ( scalar split( // , $ref ) == 1 ) {
			if ( $allele1 eq $ref ) {
				$refCount += 1;
			} 
			if ( $allele2 eq $ref ) {
				$refCount += 1;
			}
			if ( $refCount == 2 ) {
				$rstatus = "homozygous";
			} elsif ( $refCount == 1 ) {
				$rstatus = "heterozygous";
			}
		}
		$astatus = "missing";
		if ( scalar split( // , $alt ) == 1 ) {
			if ( $allele1 eq $alt ) {
				$altCount += 1;
			} 
			if ( $allele2 eq $alt ) {
				$altCount += 1;
			}
			if ( $altCount == 2 ) {
				$astatus = "homozygous";
			} elsif ( $altCount == 1 ) {
				$astatus = "heterozygous";
			}
		}
		$OUT->print( join( "\t" , ( $rstatus , $astatus , $rsIDs{$gen} , $chr , $pos , $allele1 , $allele2 , $line ) )."\n" );
	}
}
$IN2->close();
$OUT->close();

__DATA__
0	chrom
1	pos
2	ref
3	alt
4	measureset_type
5	measureset_id
6	rcv
7	allele_id
8	symbol
9	hgvs_c
10	hgvs_p
11	molecular_consequence
12	clinical_significance
13	pathogenic
14	benign
15	conflicted
16	review_status
17	gold_stars
18	all_submitters
19	all_traits
20	all_pmids
21	inheritance_modes
22	age_of_onset
23	prevalence
24	disease_mechanism
25	origin
26	xrefs
