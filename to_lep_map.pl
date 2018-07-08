#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Std;
my %opts;

##v2
unless (getopts("g:b:m:p:f:n:h", \%opts))
{
        printhelp();
        print "Error: some options are not set properly!\n";
        exit;
}

if (defined $opts{"h"}) 
{
	 printhelp();
	 exit;
}

my $familyName = "F";
if (defined $opts{"n"}) 
{
	 $familyName = $opts{"n"};
}

my $inputfile = "";
if ($opts{"g"}) 
{
	$inputfile = $opts{"g"};
}
else
{
	printhelp();
    print "Error: input genotype file is not specified!\n";
    exit;
}


my $maf = 0;
if ($opts{"f"}) 
{
	$maf = $opts{"f"};
}


my %moms = ();
my $hasmom = -1;
if ($opts{"m"}) 
{
	my @t = $opts{"m"}=~/(\d+)/g;
	$hasmom =$t[0] -1;
	foreach  (@t) 
	{
		$moms{$_ - 1} ="";
	}
}

my %dads = ();
my $hasdad = -1;
if ($opts{"p"}) 
{
	my @t = $opts{"p"}=~/(\d+)/g;
	$hasdad =$t[0] - 1;
	foreach  (@t) 
	{
		$dads{$_ - 1} ="";
	}
}

my %blanks = ();
if ($opts{"b"}) 
{
	my @t = $opts{"b"}=~/(\d+)/g;
	foreach  (@t) 
	{
		$blanks{$_ - 1} ="";
	}
}


open OUT, ">$familyName.lepmap3.pedigree.txt";
open OUT2, ">$familyName.vcf";
open OUT3, ">$familyName.vcf.allelelookup.txt";
open OUT4, ">$familyName.lepmap2.txt";
open OUT5, ">$familyName.lepmap2.markers.txt";
open OUT6, ">$familyName.rqtl.txt";


	print OUT2 "##fileformat=VCFv4.0\n";
	print OUT2"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
	print OUT2"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth.\">\n";
	#print OUT2"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
	print OUT2"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	#print OUT2"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n";
	print OUT2 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";


	print OUT3 "Locus\tAllele\tAlleleIndex\tHaplotypeId\n";

my @rqtl;
my @lepmap2;
my @lepmap2Markers;
open (IN, "$inputfile")  || die "Error: cannot open input genoytpe file $inputfile\n";

my $header = <IN>;
$header=~s/\s+$//;
my @samples = split "\t", $header;
shift @samples;
shift @samples;

my @Samples = ();
my @usedSamples = (); 
my @usedSamplesDad = (); 
my @usedSamplesMom = (); 
my @usedSamplesSex = (); 
my @usedSampleIndex = ();
my $usedSampleCount = 0;


my $momSample = 0;
my $dadSample = 0;
if ($hasdad>=0)
{
	$usedSampleCount++;
	$dadSample = $samples[$hasdad];
	unless ($dadSample) 
	{
		my $c = @samples + 0;
		print "Error: Father sample $hasdad does not exist! Total sample count is $c \n"; 
		exit;
	}
	push @usedSamples, $dadSample;
	push @usedSamplesDad, 0;
	push @usedSamplesMom, 0;
	push @usedSamplesSex, 1;
	push @usedSampleIndex, $hasdad;
}

if ($hasmom>=0)
{
	$usedSampleCount++;
	$momSample = $samples[$hasmom];
	unless ($dadSample) 
	{
		my $c = @samples + 0;
		print "Error: Mother sample $hasmom does not exist! Total sample count is $c \n"; 
		exit;
	}
	push @usedSamples, $momSample;
	push @usedSamplesDad, 0;
	push @usedSamplesMom, 0;
	push @usedSamplesSex, 2;
	push @usedSampleIndex, $hasmom;
}


LOOP3:for (my $i=0; $i<=$#samples; $i++) 
{
	if (exists $dads{$i}) 
	{
	}
	elsif (exists $moms{$i})  
	{
	}
	elsif (exists $blanks{$i})
	{
	}
	else
	{
		$usedSampleCount ++;
		push @usedSamples, $samples[$i];
		push @usedSamplesMom, $momSample;
		push @usedSamplesDad, $dadSample;
		push @usedSamplesSex, 0;
		push @usedSampleIndex, $i;
	}
}


## first line
print OUT "Contig\tPos";
for (my $i=0; $i< $usedSampleCount; $i++) 
{
	print OUT "\t".$familyName;
}


#second line
print OUT "Contig\tPos\t";
print OUT (join "\t", @usedSamples);
print OUT "\n";

#third line
print OUT "Contig\tPos\t";
print OUT (join "\t", @usedSamplesDad);
print OUT "\n";

#fourth line
print OUT "Contig\tPos\t";
print OUT (join "\t", @usedSamplesMom);
print OUT "\n";

#fifth line
print OUT "Contig\tPos\t";
print OUT (join "\t", @usedSamplesSex);
print OUT "\n";


## six line
print OUT "Contig\tPos";
for (my $i=0; $i< $usedSampleCount; $i++) 
{
	print OUT "\t0";
}
print OUT "\n";



print OUT2 (join "\t", @usedSamples);
print OUT2 "\n";


for (my $i=0; $i< $usedSampleCount; $i++) 
{
	push @{$lepmap2[$i]}, $familyName;
	push @{$lepmap2[$i]}, $usedSamples[$i];
	push @{$lepmap2[$i]}, $usedSamplesDad[$i];
	push @{$lepmap2[$i]}, $usedSamplesMom[$i];
	push @{$lepmap2[$i]}, $usedSamplesSex[$i];
	push @{$lepmap2[$i]}, 1;

	push @{$rqtl[$i]}, $usedSamples[$i];
}

my %markernameToContig;
my %markernameToPos;
LOOP1:while (<IN>) 
{
	s/\s+$//;
	my @data = split "\t";
	my $locus = shift @data;
	my $alleles =shift @data;
	my ($contig, $pos);
	if ($locus=~/(.+)_(\d+)$/) 
	{
		$contig = $1;
		$pos = $2;
	}
	else
	{
		print "warning: marker $locus is not in right format! output as contig $locus position 0bp\n";
		$contig = $locus;
		$pos = 0;
	}
	$markernameToContig{$locus} = $contig;
	$markernameToPos{$locus}= $pos;
	unless ($pos) 
	{
		$pos = "";
	}
	$alleles =~s/;\s*$//;
	my @alleles = split ";", $alleles;
	my %usedHapToVCFAlleleIndex = ();
	my @usedAllelesSorted = ();
	my @usedAllelesSortedNT = ();
	my $alleleCount =0;
	my @aindex = (0,1,2,3);
	my @aNT = qw(A C G T);
	my %aindexToNt = (0=>"A", 1=>"C", 2=>"G", 3=>"T");
	LOOP2:foreach  (@alleles) 
	{
		last LOOP2 if ($alleleCount>=4);
		if(/(\d+)\(([0-9.]+)\)/)
		{
			my ($a, $f) = ($1, $2);
			if ($f>$maf) 
			{
				my $tindex = shift @aindex; 
				my $t = shift @aNT;
				$usedHapToVCFAlleleIndex{$a} = $tindex;
				push @usedAllelesSorted, $a;
				push @usedAllelesSortedNT, $t;
				$alleleCount ++;
				print OUT3 $locus, "\t$t\t$tindex\t$a\n", 
			}

		}		
	}
	next LOOP1 if ($alleleCount<2);
	my $ref = $usedAllelesSortedNT[0];
	my $alt = $usedAllelesSortedNT[1];
	if (@usedAllelesSorted >2) 
	{
		$alt = join ",", @usedAllelesSortedNT[1..$#usedAllelesSortedNT];
	}
	
	push @lepmap2Markers, $locus;
	print OUT2 $contig, "\t", $pos, "\t", $locus, "\t", $ref, "\t", $alt, "\t.\tPASS\t.\tGT:AD:DP" ;

	my $sampleIndexInLM2 = 0;
	foreach  my $i (@usedSampleIndex) 
	{
		my ($g, $d)  = split ":",  $data[$i];
		
		my @gts = split /\//, $g;
		my @ds = split ",", $d;
		my @outgts = ();
		my @outgts_lm2 = ();
		my %alleleToDps = ();
		foreach my $gt (@gts) 
		{
			if (@ds>0) 
			{
				my $ad = shift @ds;
				if (exists $usedHapToVCFAlleleIndex{$gt}) 
				{
					my $t = $usedHapToVCFAlleleIndex{$gt}; 
					push @outgts, $t;
					push  @outgts_lm2, $gt;
					$alleleToDps{$gt} = $ad;
				}
			}

		}
		my $outgt = "./.";
		my $outgt_lm2 = "0 0";
		my $outgt_rqtl = "NA";
		if (@outgts == 2) 
		{
			@outgts = sort {$a<=>$b} @outgts;
			@outgts_lm2 = sort {$a<=>$b} @outgts_lm2;
			$outgt = $outgts[0]. "/". $outgts[1];
			$outgt_lm2 = $outgts_lm2[0]. " ". $outgts_lm2[1];
			$outgt_rqtl = $aindexToNt{$outgts[0]} . $aindexToNt{$outgts[1]};
		}
		elsif (@outgts == 1) 
		{
			$outgt = $outgts[0]. "/". $outgts[0];
			$outgt_lm2 =  $outgts_lm2[0]. " ". $outgts_lm2[0];
			$outgt_rqtl = $aindexToNt{$outgts[0]} . $aindexToNt{$outgts[0]};
		}
		
		my $ADstr = "";
		my $tDP = 0;
		foreach  my $a (@usedAllelesSorted) 
		{
			if (exists $alleleToDps{$a}) 
			{
				$ADstr .= $alleleToDps{$a}. ",";
				$tDP += $alleleToDps{$a};
			}
			else
			{																
				$ADstr .= "0,"
			}
		}
		$ADstr =~ s/,$//;
		print OUT2 "\t$outgt:$ADstr:$tDP";
		push  @{$lepmap2[$sampleIndexInLM2]}, $outgt_lm2;
		push @{$rqtl[$sampleIndexInLM2]}, $outgt_rqtl;
		$sampleIndexInLM2 ++;
		#print join ";",@usedAllelesSorted;
		#print "\t$i\t$data[$i]\t". (values %alleleToDps). "\t$outgt:$ADstr:$tDP\n";

	}
	#exit;
	print OUT2 "\n";
	 


}
close OUT;
close OUT2;
close OUT3;


foreach  ( @lepmap2Markers) 
{
	print OUT6 "\t$_";
}
print OUT6 "\n";

foreach  ( @lepmap2Markers) 
{
	print OUT6 "\t$markernameToContig{$_}";
}
print OUT6 "\n";

foreach  ( @lepmap2Markers) 
{
	print OUT6 "\t$markernameToPos{$_}";
}
print OUT6 "\n";


for (my $i=0; $i< $usedSampleCount; $i++) 
{
	print OUT4 (join "\t", @{$lepmap2[$i]});
	print OUT4 "\n";

	print OUT6 (join "\t", @{$rqtl[$i]});
	print OUT6 "\n";
}

close OUT4;



my $i=0;
foreach  (@lepmap2Markers) 
{
	print OUT5  $i, "\t", $_, "\n";
	$i ++;
}

close OUT5;


sub printhelp
{
	print "Usage: analyze_amplicon.pl -s mySampleFile -k myKeyFile\n";
	print "Options:\n";
	print "-g: Input genotype file\n";
	print "-f: mininum allele frequency\n";
	print "-b: define blank samples index, if there are several separated by comma. First sample index is 1.\n";
	print "-m: mather parent index, if there are several separated by comma. First sample index is 1.\n";
	print "-p: paternal parent index, if there are several separated by comma. First sample is 1.\n";
	print "-n: family name";
	print "-h: help menu.";
}