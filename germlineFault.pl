#!/usr/bin/perl
#11 February 2017 - Adam D Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

my $usage = 'perl germlineFault.pl <child_genome> <mother_genome> <father_genome> <child_gender> <output> 
';

die $usage , unless @ARGV == 5;
my ( $childGenome , $motherGenome , $fatherGenome , $gender , $output ) = @ARGV;

my $IN1 = FileHandle->new( "$childGenome" , "r" );
if ( not defined $IN1 ) { die "ADSERROR: Could not open/read $childGenome\n"; }

my $IN2 = FileHandle->new( "$motherGenome" , "r" );
if ( not defined $IN2 ) { die "ADSERROR: Could not open/read $motherGenome\n"; }

my $IN3 = FileHandle->new( "$fatherGenome" , "r" );
if ( not defined $IN3 ) { die "ADSERROR: Could not open/read $fatherGenome\n"; }

my $OUT = FileHandle->new( "$output" , "w" );
if ( not defined $OUT ) { die "ADSERROR: Could not open/write $output\n"; }

my $snps;
my %rsIDs;
while ( my $line = $IN1->getline ) {
	next if ( $line =~ m/^#/ );
	chomp( $line );
	my ( $rsID , $chr , $pos , $genotype ) = split( "\t" , $line );
	my ( $allele1 , $allele2 ) = split // , $genotype;
	if ( $chr eq "MT" ) { $allele2 = "."; }
	if ( $gender eq "m" ) {
		if ( $chr eq "X" or $chr eq "Y" ) {
			$allele2 = ".";
		}
	}
	my $gen = join( ":" , ( $chr , $pos ) );
 	$snps->{$gen}->{$allele1}->{"mother"} = -1;
 	$snps->{$gen}->{$allele2}->{"mother"} = -1;
 	$snps->{$gen}->{$allele1}->{"father"} = -1;
 	$snps->{$gen}->{$allele2}->{"father"} = -1;
	$rsIDs{$gen} = $rsID;
}
$IN1->close();

&matchParent( $IN2 , $snps , "mother" );
&matchParent( $IN3 , $snps , "father" );

$OUT->print( "dbSNP\tGenomicPosition\tAlleles\tFromMother\tFromFather\n" );
foreach my $gen ( sort keys %{$snps} ) {
	my @fromMother;
	my @fromFather;
	my @alleles;
	foreach my $allele ( keys %{$snps->{$gen}} ) {
		push @alleles , $allele;
		if ( $snps->{$gen}->{$allele}->{"mother"} == -1 ) {
			push @fromMother , "-";
		} else {
			push @fromMother , $snps->{$gen}->{$allele}->{"mother"};
		}
		if ( $snps->{$gen}->{$allele}->{"father"} == -1 ) {
			push @fromFather , "-";
		} else {
			push @fromFather , $snps->{$gen}->{$allele}->{"father"};
		}
	}
	$OUT->print( join( "\t" , ( $rsIDs{$gen} , $gen , join( "," , @alleles ) , join( ":" , @fromMother ) , join( ":" , @fromFather ) ) )."\n" );
}
	

sub matchParent {
	my ( $fh , $snps , $parent ) = @_;
	while ( my $line = $fh->getline ) {
		next if ( $line =~ m/^#/ );
		chomp( $line );
		my ( $rsID , $chr , $pos , $genotype ) = split( "\t" , $line );
		my ( $allele1 , $allele2 ) = split // , $genotype;
		if ( $chr eq "MT" ) { $allele2 = "."; }
		my $gen = join( ":" , ( $chr , $pos ) );
		if ( exists $snps->{$gen} ) {
			if ( exists $snps->{$gen}->{$allele1} ) {
				$snps->{$gen}->{$allele1}->{$parent} = 1;
			}
			if ( exists $snps->{$gen}->{$allele2} ) {
				$snps->{$gen}->{$allele2}->{$parent} = 1;
			}
		}
	}
	$fh->close();
}
