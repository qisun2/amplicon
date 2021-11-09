#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
my %opts;

##v2
unless (getopts("g:b:m:p:j:k:f:n:l:d:xh", \%opts))
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
my $forcetop4 = 0;
if (defined $opts{"x"}) 
{
	 $forcetop4 = 1;
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

my $mindepth = 0;
if ($opts{"d"}) 
{
	$mindepth = $opts{"d"};
}


my %moms = ();
my @moms = ();
my $hasmom = -1;
my $mname = "";
if ($opts{"m"}) 
{
	my @t = $opts{"m"}=~/(\d+)/g;
	foreach  (@t) 
	{
		$moms{$_ - 1} ="";
	}
	@moms = sort {$a<=>$b} keys %moms;
	$hasmom =$moms[0];

	if ($opts{"j"}) 
	{
		$mname = $opts{"j"};
	}

}

my %dads = ();
my @dads = ();
my $hasdad = -1;
my $pname = "";
if ($opts{"p"}) 
{
	my @t = $opts{"p"}=~/(\d+)/g;
	foreach  (@t) 
	{
		$dads{$_ - 1} ="";
	}
	@dads = sort {$a<=>$b} keys %dads;
	$hasdad =$dads[0];
	if ($opts{"k"}) 
	{
		$pname = $opts{"k"};
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


my %marker2chrpos ;
my $lookuptable = 0;
if ($opts{"l"}) 
{
	$lookuptable = 1;
	my $markerlookupfile  = $opts{"l"};
	open (IN, "$markerlookupfile")  || die "Error: cannot open marker to chromosome position file $markerlookupfile\n";
	while (<IN>) 
	{
		chomp;
		my ($marker, $chr, $pos) = split "\t";
		$marker2chrpos{$marker} = "$chr\t$pos";
	}
	close IN;
	
}


##if both mom and dad available use alleles of mom and dad
my $usealleles = "top4";
if (($hasmom>=0) && ($hasdad>=0) && ($forcetop4==0)) 
{
	$usealleles = "parents";
}


open OUT, ">$familyName.lepmap3.pedigree.txt";
open OUT2, ">$familyName.vcf";
open OUT3, ">$familyName.vcf.allelelookup.txt";

#open OUT4, ">$familyName.lepmap2.txt";
#open OUT5, ">$familyName.lepmap2.markers.txt";
#open OUT6, ">$familyName.rqtl.txt";


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
	if ($pname ne "") 
	{
		$dadSample  = $pname;
	}
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

	if ($mname ne "") 
	{
		$momSample  = $mname;
	}

	unless ($momSample) 
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


### first line
print OUT "Contig\tPos";
for (my $i=0; $i< $usedSampleCount; $i++) 
{
	print OUT "\t".$familyName;
}
print OUT "\n";

##second line
print OUT "Contig\tPos\t";
print OUT (join "\t", @usedSamples);
print OUT "\n";
#
##third line
print OUT "Contig\tPos\t";
print OUT (join "\t", @usedSamplesDad);
print OUT "\n";
#
##fourth line
print OUT "Contig\tPos\t";
print OUT (join "\t", @usedSamplesMom);
print OUT "\n";
#
##fifth line
print OUT "Contig\tPos\t";
print OUT (join "\t", @usedSamplesSex);
print OUT "\n";
#
#
### six line
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
my %outputLines = ();
my %outputLinesPos = ();

my $unknownPosMarker_index = 0;
LOOP1:while (<IN>) 
{
	s/\s+$//;
	my @data = split "\t";
	my $locus = shift @data;
	my $alleles =shift @data;
	my ($contig, $pos);
	if ($lookuptable) 
	{
		unless (exists $marker2chrpos{$locus}) 
		{
			next LOOP1;
		}
		($contig, $pos) = split "\t", $marker2chrpos{$locus};
	}
	elsif ($locus=~/(.+)_(\d+)$/) 
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
		$contig = "unknown";
		$unknownPosMarker_index  = $unknownPosMarker_index  + 100;
		$pos = $unknownPosMarker_index ;
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


	## merging

	merge_indv(\@moms, \@data);


	merge_indv(\@dads, \@data);



	if ($usealleles eq "top4") 
	{
		LOOP2:foreach  (@alleles) 
		{
			last LOOP2 if ($alleleCount>=4);
			if(/(\d+)\(([0-9.]+)\)/)
			{
				my ($a, $f) = ($1, $2);
				#print "fff$f $maf\n";
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
	}
	elsif ($usealleles eq "parents")
	{
		my @parents = ($hasmom, $hasdad);
		my %checkParentAlleles = ();

		LOOPP:foreach my $i (@parents) 
		{
			my ($g, $d)  = split ":",  $data[$i];
			next LOOPP unless ($g=~/\d/);
			my @gts = split /\//, $g;
			my @ds = split ",", $d;

			if ($gts[0] == $gts[1]) 
			{
				$checkParentAlleles{$gts[0]} += $ds[0];
			}
			else
			{
				$checkParentAlleles{$gts[0]} += $ds[0];
				$checkParentAlleles{$gts[1]} += $ds[1];
			}

		}
		#$,="\t";
		#print %checkParentAlleles, "\n";
		
		@alleles = reverse sort {$checkParentAlleles{$a}<=>$checkParentAlleles{$b}} (keys %checkParentAlleles);

		LOOPP2:foreach my $a  (@alleles) 
		{
			last LOOPP2 if ($alleleCount>=4);
			my $tindex = shift @aindex; 
			my $t = shift @aNT;
			$usedHapToVCFAlleleIndex{$a} = $tindex;
			push @usedAllelesSorted, $a;
			push @usedAllelesSortedNT, $t;
			$alleleCount ++;
			print OUT3 $locus, "\t$t\t$tindex\t$a\n", 
	
		}
	}

	my ($ref, $alt);
	if ($alleleCount==1)
	{
		$ref = $usedAllelesSortedNT[0];
	    $alt = '.';
	}
	elsif ($alleleCount >1)
	{
		$ref = $usedAllelesSortedNT[0];
	    $alt = $usedAllelesSortedNT[1];
	}
	else
	{
		next LOOP1;
	}

	if (@usedAllelesSorted >2) 
	{
		$alt = join ",", @usedAllelesSortedNT[1..$#usedAllelesSortedNT];
	}
	
	push @lepmap2Markers, $locus;
	my $outline = join "\t", ($contig, $pos, $locus, $ref,  $alt, ".\tPASS\t.\tGT:AD:DP" );
	#print OUT2 $contig, "\t", $pos, "\t", $locus, "\t", $ref, "\t", $alt, "\t.\tPASS\t.\tGT:AD:DP" ;

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

		if ($tDP<$mindepth) {
			$outgt = "./.";
		}
		$outline .= "\t$outgt:$ADstr:$tDP";
		#print OUT2 "\t$outgt:$ADstr:$tDP";
		push  @{$lepmap2[$sampleIndexInLM2]}, $outgt_lm2;
		push @{$rqtl[$sampleIndexInLM2]}, $outgt_rqtl;
		$sampleIndexInLM2 ++;
		#print join ";",@usedAllelesSorted;
		#print "\t$i\t$data[$i]\t". (values %alleleToDps). "\t$outgt:$ADstr:$tDP\n";

	}
	#exit;
	push @{$outputLines{$contig}}, $outline;
	push @{$outputLinesPos{$contig}}, $pos;


	#print OUT2 "\n";
}

my @contigs = sort keys %outputLines;
foreach my $contig (@contigs) 
{
	my @aa = @{$outputLinesPos{$contig}};
	my @idx = sort { $aa[$a] <=> $aa[$b] } 0 .. $#aa;
	foreach  (@idx) 
	{
		print OUT2 $outputLines{$contig}[$_];
		print OUT2 "\n"
	}
}

#close OUT;
close OUT2;
close OUT3;


## sort by genotyping results


#foreach  ( @lepmap2Markers) 
#{
#	print OUT6 "\t$_";
#}
#print OUT6 "\n";
#
#foreach  ( @lepmap2Markers) 
#{
#	print OUT6 "\t$markernameToContig{$_}";
#}
#print OUT6 "\n";
#
#foreach  ( @lepmap2Markers) 
#{
#	print OUT6 "\t$markernameToPos{$_}";
#}
#print OUT6 "\n";
#
#
#for (my $i=0; $i< $usedSampleCount; $i++) 
#{
#	print OUT4 (join "\t", @{$lepmap2[$i]});
#	print OUT4 "\n";
#
#	print OUT6 (join "\t", @{$rqtl[$i]});
#	print OUT6 "\n";
#}
#
#close OUT4;



#my $i=0;
#foreach  (@lepmap2Markers) 
#{
#	print OUT5  $i, "\t", $_, "\n";
#	$i ++;
#}
#
#close OUT5;

sub merge_indv
{
	my $merginglist = shift;
	my @merginglist = @{$merginglist};
	my $firstmember = $merginglist[0];

	my $gtdata = shift;
	my @gtdata = @{$gtdata};

	my %gt2counts = ();
	LOOPTTT:foreach my $i (@merginglist) 
	{
		my ($g, $d)  = split ":",  $gtdata[$i];
		if ($d eq "0") 
		{
			next LOOPTTT;
		}
		
		my @gts = split /\//, $g;
		my @ds = split ",", $d;
		if (@ds ==1) 
		{
			$gt2counts{$gts[0]} += $ds[0];
		}
		else
		{
			foreach my $gt (@gts) 
			{
				my $count = shift @ds;
				$gt2counts{$gt} += $count;
			} 
		}

	}
	#$,="\t";
	#print %gt2counts, "\n";
	my @mergedgts = reverse sort {$gt2counts{$a}<=>$gt2counts{$b}} keys %gt2counts;
	if (@mergedgts == 0) 
	{
		return;
	}
	elsif (@mergedgts ==1) 
	{
		my $a = $mergedgts[0];
		my $c = $gt2counts{$a};
		${$gtdata}[$firstmember] = "$a/$a:$c";
	}
	else
	{
		my $a = $mergedgts[0];
		my $b = $mergedgts[1];
		my $ac = $gt2counts{$a};
		my $bc = $gt2counts{$b};
		${$gtdata}[$firstmember] = "$a/$b:$ac,$bc";
	}
}

sub printhelp
{
	print "Usage: to_lep_map.pl -g inputGTFile -f maf -b <blank samples> -m  <maternal samples> -p <paternal samples> -n familyName\n";
	print "Options:\n";
	print "-g: Input genotype file, the hap_genotype file from amplicon.py\n";
	print "-f: mininum allele frequency\n";
	print "-b: define skipped samples index, if there are several separated by comma. First sample index is 1.\n";
	print "-m: maternal parent index, if there are several separated by comma. First sample index is 1.\n";
	print "-p: paternal parent index, if there are several separated by comma. First sample is 1.\n";
	print "-j: maternal parent name, if not set, use the first marternal individual name.\n";
	print "-k: paternal parent name, if not set, use the first paternal individual name.\n";
	print "-n: family name\n";
	print "-l: provide a table to convert haplotype marker to chromosome name and position. A tab-delimited table with 3 columns: hapmarker, chr, pos\n";
	print "-x: force using top 4 alleles even if parents are present\n";
	print "-d: minimum genotype read depth. genotype below this depth will be converted to unknown\n";
	print "-h: help menu.\n";
}
