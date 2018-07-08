#!/usr/bin/perl
use strict;
use warnings;
use Spreadsheet::ParseXLSX;
use Excel::Writer::XLSX;
use Getopt::Std;

my %opts;
my $excelinput;
my $allelefile;
my $outputfile;
my $outdir;
my $mode;
my $tmpdir = "/tmp";
my $max_mismatch = 0;
our $fragsize = 100;
## flash parameter
my $threads = 8;
my $min_overlap_length = 20;
my  $max_overlap_length = 350;
my $contiging_mismatch = 0.03;
our $flash_cmd = "/programs/FLASH2/flash2 -t $threads -q -m $min_overlap_length -M $max_overlap_length -x $contiging_mismatch ";

our $bwacmd = "bwa mem -t $threads -M -v 1 ";


our %compnt = ("A"=>"T",
                        "T"=>"A",
                        "G"=>"C",
                        "C"=>"G",
                        "N"=>"N",
                        );


##verify input parameters
unless (getopts("i:o:m:r:a:t:x:s:h", \%opts))
{
        & print_help();
        print "Error: some options are not set properly!\n";
        exit;
}

if (defined $opts{"h"}) 
{
	 & print_help();
	 exit;
}

if (defined $opts{"x"}) 
{
	 $max_mismatch = $opts{"x"};
	 unless ($max_mismatch =~/^\d+$/) 
	 {
		& print_help();
        print "Error: -x max mismatch must be an integer!\n";
        exit;
	 }
}

if (defined $opts{"t"}) 
{
	 $tmpdir = $opts{"t"};
	 unless (-d $tmpdir) 
	 {
		 mkdir $tmpdir;
	 }
}

our $cutadapt_cmd = "cutadapt --quiet -a RRR --trimmed-only XXX  | cutadapt --quiet -g FFF --trimmed-only - > $tmpdir/trimmed.fastq\n";


if (defined $opts{"i"}) 
{
	$excelinput = $opts{"i"};
}
else
{
	& print_help();
    print "Error: input excel file is required!\n";
    exit;
}


if (defined $opts{"o"}) 
{
	$outputfile = $opts{"o"};
}
else
{
	& print_help();
    print "Error: output excel file name is required!\n";
    exit;
}

if (defined $opts{"r"}) 
{
	$outdir = $opts{"r"};
	unless (-d $outdir) 
	{
		mkdir $outdir;
	}
}


if (defined $opts{"s"}) 
{
	$fragsize = $opts{"s"};
	unless ($fragsize=~/^\d+$/) 
	{
		& print_help();
		print "Error: minimum size must be a number!\n";
		exit;

	}
}

if (defined $opts{"a"}) 
{
	$allelefile = $opts{"a"};
	unless (-e $allelefile) 
	{
		print "Error: Cannot find the allele file $allelefile.\n";
		exit;
	}
}


if (defined $opts{"m"}) 
{
	$mode = $opts{"m"};
}
else
{
	& print_help();
    print "Error: you must specify running mode! The acceptable values are '-m d' or '-m p'\n";
    exit;
}

if ($mode=~/^d/i) 
{
	unless (defined $outdir) 
	{
		& print_help();
		print "Error: you use '-m d' for discovery mode with no known alleles. You must specify an output directory for this mode use '-r directoryName'\n";
		exit;
	}
}
elsif ($mode=~/^c/i) 
{
	unless (defined $allelefile) 
	{
		& print_help();
		print "Error: you use '-m c' for counting known allelels. You must specify a haplotype allele files with  '-a hapAlleleFile'\n";
		exit;
	}
}
elsif ($mode=~/^p/i)
{
}
else
{
	& print_help();
    print "Error: you must specify running mode! The acceptable values are 'd' 'p' 'c'.\n";
    exit;
}


unless ($outputfile=~/xlsx$/) 
{
	& print_help();
	print "Error: The output file name you entered is $outputfile. The file name must end with '.xlsx'!";
	exit;
}

unless ($excelinput =~/xlsx$/) 
{
	& print_help();
	print "Error: The input file name you entered is $excelinput. An excel spreadsheet *.xlsx file is required!";
	exit;
}

unless (-e $excelinput) 
{
	print "Error: Cannot find the file $excelinput.\n";
	exit;
}






# quality filter
our $minQ = 10;
$minQ += 33; ## ascii value for fastq quality 



unless (-e $excelinput) 
{
	print "Cannot find the file $excelinput.\n";
	exit;
}

unless (-d $tmpdir) 
{
	mkdir $tmpdir;
}



if (-e "$outputfile") 
{
	print "There is an existing file $outputfile.";
	print "Please delete or re-naming the existing file first\n";
	exit;
}

## read input excel file
my $parser   = Spreadsheet::ParseXLSX->new();
my $workbook = $parser->parse($excelinput);
my @worksheets = $workbook->worksheets();

our %trait2marker = ();
our @traits = ();

our %family2ind = ();

our %family2traits = ();
our @families = ();

our %ind2file = ();

our @inds;
our @parents1;
our @parents2;
our %indtype;

our @familytraits = ();
our @alltraits = ();
our @alldonors = ();
our @allmarkers = ();
our @allprimer1 = ();
our @allprimer2 = ();
our @all_d_alleles = ();
our @all_u_alleles = ();
our $total_markers = 0;
our %marker2primers = ();


print "Step 1. Read spreadsheet.\n";
foreach my $sheet (@worksheets) 
{
	&read_sheet ($sheet);
}

print "Step 2. Evaluate spreadsheet data.\n";
unless (@traits>0) 
{
	print "Error: either breeding_markers sheet is missing or breeding_markers sheet has no values\n";
	exit;
}



unless (((keys %family2ind)+0)>0) 
{
	print "Error: either family sheet is missing or family sheet has no values\n";
	print keys %family2ind;
	exit;
}

unless (@families>0) 
{
	print "Error: either family_trait sheet is missing or family_trait sheet has no values\n";
	exit;
}

unless ((keys %ind2file)>0) 
{
	print "Error: either Seq_file sheet is missing or family Seq_file has no values\n";
	exit;
}


## evaluate familiies
foreach my $FamilyName (@families) 
{
	my @familytraits = @{$family2traits{$FamilyName}};
	foreach  (@familytraits) 
	{
		my ($t, $d) = split "\t";
		unless (exists $trait2marker{$t}) 
		{
			print "Family $FamilyName has a trait $t, but this trait $t does not have markers\n";
			exit;
		}
	}

	unless (exists $family2ind{"progeny"}{$FamilyName}) 
	{
		print "Family $FamilyName has no progenies.\n";
		exit;
	}

	@inds = @{$family2ind{"progeny"}{$FamilyName}};
	@parents1 = @{$family2ind{"parent1"}{$FamilyName}};
	@parents2 = @{$family2ind{"parent2"}{$FamilyName}};



	@inds = ( @parents1,  @parents2, @inds);
	foreach my $ind (@inds) 
	{
		unless (exists $ind2file{$ind}) 
		{
			print "Warning: individual ${ind} in family $FamilyName does not have file name\n";
		}
	}
}

foreach my $trait (@traits) 
{
	if (exists $trait2marker{$trait}) 
	{
		my @markers = @{$trait2marker{$trait}};
		foreach  (@markers) 
		{
			my ($markerid, $primer1, $primer2, $dalleles, $ualleles)=split "\t";
			push @allmarkers, $markerid;
		}
	}
	else
	{
		print "Error: trait $trait  does not have markers in the breeding_markers sheet of the input excel file\n";
	    exit;
	}
}



## evaluate allele file and create allele database
our %alleleHash=();
my %markerInAlleleFile = ();
if (defined $allelefile) 
{
	open (IN, $allelefile) || die "cannot open $allelefile";

	while (<IN>) 
	{
		if (/\w/) 
		{
			s/\s+$//;
			my ($marker, $id, $seq) = split "\t";
			$markerInAlleleFile{$marker} = "";
			unless ($seq=~/\w/) 
			{
				print "Error: in $allelefile, marker $marker has no sequence\n";

			}
			$seq=~s/\W//g;

			$alleleHash{$seq} = "${marker}\t${id}";
		}
	}
	close IN;

	if ($mode=~/^c/i) 
	{
		my $e = 0;
		my $u =0;
		my @u = ();
		foreach  (@allmarkers) 
		{
			if (exists $markerInAlleleFile{$_}) 
			{
				$e ++;
			}
			else
			{
				$u ++;
				push @u, $_;
			}
		}
		print "For all markes to be scored, $e exist in $allelefile; $u does not.";
		if ($u>0) 
		{
			print "Please check markers in the breeding_markers spreadsheet, and modify/delete these markers: \n";
			print join "\n", @u;
			print "\n";
			exit;
		}
	}

}

print "Step 3. Process families.\n";
our $newworkbook  = Excel::Writer::XLSX->new( $outputfile);
our $newsheet;
our @rowcells; 
our $rowcount;

our $green_gb = $newworkbook->add_format();
our $mom_gb = $newworkbook->add_format();
our $dad_gb = $newworkbook->add_format();
our $regular_gb = $newworkbook->add_format();
our $bold = $newworkbook->add_format();

$bold->set_bold(); 

$green_gb->set_bg_color( "#8fbc8f" );
$mom_gb->set_bg_color( "#eaadea" );
$dad_gb->set_bg_color( "#99ccff" );





## process families
our $FamilyName;
FAMILYLOOP:foreach $FamilyName (@families) 
{
	print "Families: $FamilyName\n";
	unless (exists $family2traits{$FamilyName}) 
	{
		print "Family $FamilyName has no traits to score.\n";
		next FAMILYLOOP;
	}

	@familytraits = @{$family2traits{$FamilyName}};
	@alltraits = ();
	@alldonors = ();
	@allmarkers = ();
	@allprimer1 = ();
	@allprimer2 = ();
	@all_d_alleles = ();
	@all_u_alleles = ();
	$total_markers = 0;
	foreach  (@familytraits) 
	{
		my ($trait, $donor) = split "\t";
		my @markers = @{$trait2marker{$trait}};
		foreach  (@markers) 
		{
			my ($markerid, $primer1, $primer2, $dalleles, $ualleles)=split "\t";
			push @alltraits, $trait;
			push @alldonors, $donor;
			push @allmarkers, $markerid;
			push @allprimer1, $primer1;
			push @allprimer2, $primer2;
			push @all_d_alleles, $dalleles; 
			push @all_u_alleles, $ualleles;
		}
	}
	$total_markers = @allmarkers+0;
	#open FAM, ">$outputdir/$FamilyName";
	$newsheet = $newworkbook->add_worksheet( $FamilyName ); 
	@rowcells = ("", "", "Trait",  @alltraits); $newsheet->write_row(0, 0, \@rowcells);
	@rowcells = ("", "", "Donor",  @alldonors); $newsheet->write_row(1, 0, \@rowcells);
	@rowcells = ("Type", "Ind.", "Reads",  @allmarkers); $newsheet->write_row(2, 0, \@rowcells, $bold);
	$rowcount = 2;

	@inds = @{$family2ind{"progeny"}{$FamilyName}};
	@parents1 = @{$family2ind{"parent1"}{$FamilyName}};
	@parents2 = @{$family2ind{"parent2"}{$FamilyName}};

	if (@parents1>0) {print "Mother: $parents1[0]\n";}
	if (@parents2>0) {print "Father: $parents2[0]\n";}
	print "Progenies: ", @inds+0, "\n";
	%indtype = ();
	foreach  (@inds) 
	{
		$indtype{$_} = "progeny";
	}
	foreach  (@parents1) 
	{
		$indtype{$_} = "mother";
	}
	foreach  (@parents2) 
	{
		$indtype{$_} = "father";
	}


	@inds = ( @parents1,  @parents2, @inds);


	if ($mode=~/^p/i) 
	{
		&production(1);
	}
	elsif ($mode=~/^d/i) 
	{
		&discovery();
	}
	elsif ($mode=~/^c/i) 
	{
		&production(2);
	}
	else 
	{
		print "Error: mode $mode is unknown\n ";
		exit;
	}


}
$newworkbook->close();



sub discovery
{
	## prepare input files for analyze amplicon
	## sample files
	my $samplefile = "$tmpdir/sample";
	my $keyfile = "$tmpdir/key";
	my $tmpout = "$tmpdir/out";
	if (-d "$tmpdir/out") 
	{
		system ("rm -fr $tmpdir/out");
	}

	open OUT, ">$samplefile";
	my $file_counter =0;
	foreach my $ind (@inds) 
	{		
		if (exists $ind2file{$ind}) 
		{
			my @files = @{$ind2file{$ind}};
			foreach my $filestr (@files) 
			{
				my ($dir, $file1, $file2) = split "\t", $filestr;
				$file_counter ++;
				
				if ($file2=~/\w/) 
				{
					print OUT "${ind}_$file_counter\t$dir/$file1\t$dir/$file2\n";
				}
				elsif ($file1=~/\w/) 
				{
					print OUT "${ind}_$file_counter\t$dir/$file1\n";
				}
			}
		}
	}
	close OUT;

	open OUT, ">$keyfile";
	foreach my $marker (@allmarkers) 
	{
		my $primer1 = shift @allprimer1;
		my $primer2 = shift @allprimer2;
		print OUT "${marker}_F\t$primer1\n";
		if ($primer2=~/\w/) 
		{
			print OUT "${marker}_R\t$primer2\n";
		}
	}
	close OUT;

	# run analyz amplicon
	system ("analyze_amplicon.pl -i 8 -s $samplefile -k $keyfile -o $tmpout -d 20:$fragsize -v 20:350:0.05 -t $threads -m clustalo:clustal");
	system ("mkdir $outdir/$FamilyName");
	system ("cp -r $tmpout/stats $outdir/$FamilyName/ ");
	system ("cp -r $tmpout/haplotype2fasta $outdir/$FamilyName/ ");
	system ("cp $tmpout/hap_genotype $outdir/$FamilyName/ ");
	system ("cp $tmpout/hap_genotype_matrix $outdir/$FamilyName/ ");


	my %hapmatrix = ();
	my %sample2total = ();
	
	open IN, "$outdir/$FamilyName/stats";
	while (<IN>) 
	{
		chomp;
		my ($indname, $total) = split "\t";
		$sample2total{$indname} = $total;

	}
	close IN;
	open IN, "$outdir/$FamilyName/hap_genotype";
	my $header = <IN>; $header=~s/\s+$//;
	my @indInTable = split "\t", $header;
	shift @indInTable;
	shift @indInTable;

	
	
	while (<IN>) 
	{
		s/\s+$//;
		my @data = split "\t";
		my $marker = shift @data;
		my $haplotypes = shift @data;
		for (my $i=0; $i<=$#indInTable; $i++) 
		{
			$hapmatrix{$indInTable[$i]}{$marker} = $data[$i];
		}
	}
	close IN;
	foreach my $ind (@indInTable) 
	{
		my $display_ind = $ind;
		$display_ind =~s/_\d+$//;
		my $type = $indtype{$display_ind}; 
		@rowcells = ($type, $display_ind, $sample2total{$ind});
		my $gb = $regular_gb;
		if ($type eq "mother") {$gb = $mom_gb;}
		elsif ($type eq "father") {$gb = $dad_gb;}
		$rowcount ++;
		$newsheet->write_row($rowcount, 0, \@rowcells, $gb);
		for (my $i=0; $i<=$#allmarkers; $i++) 
		{
			my $marker = $allmarkers[$i];
			my $token = $hapmatrix{$ind}{$marker};
			$newsheet->write($rowcount, $i+3, $token);
		}
	}

}

sub production
{
	my $scoring_mode = shift;
	my $file_counter =0; 
	my $lasttype = "";
	foreach my $ind (@inds) 
	{		
		my $type = $indtype{$ind};
		if ($type ne $lasttype) 
		{
			$rowcount ++;
			$lasttype = $type;
		}
		print "  $type $ind: file(s)#";
		if (exists $ind2file{$ind}) 
		{
			my @files = @{$ind2file{$ind}};
			print @files+0, "\n";
			foreach my $filestr (@files) 
			{
				my ($dir, $file1, $file2) = split "\t", $filestr;
				$file_counter ++;

				print "  SampleFile #$file_counter\n";
				
				my ($total, $allelecounts) = (0, "");
				my $processfastqfile;
				if ($file2=~/\w/) 
				{
					my $flash_exec_command = "$flash_cmd -d $tmpdir -o s $dir/$file1 $dir/$file2";
					#print $flash_exec_command, "\n";
					system ("$flash_exec_command");

					$processfastqfile = "$tmpdir/s.extendedFrags.fastq";
					#($total, $allelecounts) = process_fastq("$tmpdir/s.extendedFrags.fastq", \@all_d_alleles);
				}
				elsif ($file1=~/\w/) 
				{
					$processfastqfile = "$dir/$file1";
					#($total, $allelecounts) = process_fastq("$dir/$file1", \@all_d_alleles);
				}

				if ($scoring_mode == 1) 
				{
					($total, $allelecounts) = process_fastq($processfastqfile, \@all_d_alleles);
				}
				else
				{
					($total, $allelecounts) = process_fastq_ac($processfastqfile);
					#system ("cp $tmpdir/sample.sam $ind.sam");
				}


				@rowcells = ($type, $ind, $total);
				my $gb = $regular_gb;
				if ($type eq "mother") {$gb = $mom_gb;}
				elsif ($type eq "father") {$gb = $dad_gb;}

				$newsheet->write_row($rowcount, 0, \@rowcells, $gb);

				
				for (my $i=0; $i<$total_markers;$i++) 
				{
					if ($scoring_mode == 1)
					{
						my $value = $allelecounts->[$i];
						my $normvalue = 0;
						if ($total>0) 
						{
							$normvalue = sprintf("%.0f", 1000000 * $value/$total) ;
						}
						my $token = "$value  $normvalue";
						if ($value>10) 
						{
							$newsheet->write($rowcount, $i+3, $token, $green_gb);
						}
						else
						{
							$newsheet->write($rowcount, $i+3, $token);
						}
					}
					else
					{
						my $marker = $allmarkers[$i];
						my $token = "";
						
						if (exists $allelecounts->{$marker}) 
						{
							my %allele2count = %{$allelecounts->{$marker}};
							my @alleles = reverse sort {$allele2count{$a} <=> $allele2count{$b}} (keys %allele2count);
							my $reportAlleleNum = 4;
							my $t=0;

							my $topreadcount =0;
							if (@alleles>0) 
							{
								$topreadcount = $allele2count{$alleles[0]};
							}
							ALLLOP:foreach my $allele (@alleles) 
							{
								if (($topreadcount>0) && ($allele2count{$allele}/$topreadcount <0.1)) 
								{
									last ALLLOP
								}

								$token .= "$allele:$allele2count{$allele}; ";


								$t ++;
								if ($t>=$reportAlleleNum) 
								{
									last ALLLOP;
								}
								
							}
							$token =~s/; $//;
						}
						else
						{
							$token = "./.";
						}
						$newsheet->write($rowcount, $i+3, $token);
						
					}

					
				}

				$rowcount ++ ;	

			}

		}
	}
	system ("rm $tmpdir/s.extendedFrags.fastq");
}

sub read_sheet
{
	my $sheet = shift;
	my $name = $sheet->get_name;
	$name=~s/\s//g;

	my ( $row_min, $row_max ) = $sheet->row_range();
	my ($col_min, $col_max ) = $sheet->col_range();
	print "  ${name}\n";
	
	my @data = ();
	
	for (my $row=$row_min; $row<=$row_max; $row++) 
	{
		for (my $col=$col_min; $col<=$col_max; $col++) 
		{
			my $cell = $sheet->get_cell( $row, $col );
			if (defined $cell) 
			{
				$data[$row][$col] = $cell->value();
			}
		}
	}

	if ($name=~/breeding_markers/i) 
	{
		unless (($data[0][0] =~/Trait/i) && ( $data[0][1] =~/Marker/i)  && ( $data[0][2] =~/Primer1/i) && ( $data[0][3] =~/Primer2/i) && ( $data[0][4] =~/Desirable/i) && ( $data[0][5] =~/Undesirable/i) )
		{
			print "Error: In the sheet 'breeding_markers', the top rows must have the value 'Trait' 'Marker' 'Primer1' 'Primer2' 'DesirableAlleles' 'UndesirableAlleles'\n";
			print "You have $data[0][0] $data[0][1] $data[0][2] $data[0][3] $data[0][4] $data[0][5]\n";
			exit;
		}
		for (my $row=1; $row<=$#data; $row++)
		{
			my ($trait, $marker, $primer1, $primer2, $dalleles, $ualleles) = @{$data[$row]};

			if ((defined $trait) && ($trait=~/\w/) )
			{
				unless (defined $primer1) {$primer1 = "";}
				unless (defined $primer2) {$primer2 = "";}			
				unless (defined $dalleles) {$dalleles = "";}
				unless (defined $ualleles) {$ualleles = "";}

				$marker2primers{$marker} = "$primer1\t$primer2";

				$trait=~s/\s//g; 
				$marker=~s/\s//g;
				$primer1=~s/\s//g;
				$primer2=~s/\s//g;
				$dalleles=~s/\s//g;
				$ualleles=~s/\s//g;

				if ($marker=~/\w/)
				{
					unless (exists $trait2marker{$trait}) 
					{
						push @traits, $trait;
					}
					push @{$trait2marker{$trait}}, "$marker\t$primer1\t$primer2\t$dalleles\t$ualleles";
				}
				else
				{
					print "In excel breeding_markers sheet, trait '$trait' does not have marker or DesirableAllele value\n";
					exit;
				}
			}

		}
	}
	elsif ($name=~/^family$/i)  
	{
		
		unless (( $data[0][0] =~/IndividualName/i) && ( $data[0][1] =~/Parent1/i)  && ( $data[0][2] =~/Parent2/i) && ( $data[0][3] =~/FamilyName/i))
		{
			print "Error: In the sheet 'family', the top rows must have the value 'IndividualName' 'Parent1' 'Parent2' 'FamilyName' \n";
			exit;
		}

		my %familyParent = ();
		for (my $row=1; $row<=$#data; $row++)
		{
			my ($IndividualName, $Parent1, $Parent2, $FamilyName) = @{$data[$row]};

			if ((defined $IndividualName) && ($IndividualName=~/\w/) )
			{
				$IndividualName=~s/\s//g; 
				if ($IndividualName=~/_[A-H]\d\d$/) 
				{
					$IndividualName=~s/_[A-H]\d\d$//;
				}
				$Parent1=~s/\s//g;
				$Parent2=~s/\s//g;
				$FamilyName=~s/\s//g;
				
				if ($FamilyName=~/\w/) 
				{
					unless (exists $family2ind{"parent1"}{$FamilyName}) 
					{
						if ((defined $Parent1) && ($Parent1=~/\w/)) {$Parent1=~s/_[A-H]\d\d$//; push @{$family2ind{"parent1"}{$FamilyName}}, $Parent1;}
						if ((defined $Parent2) && ($Parent2=~/\w/)) {$Parent2=~s/_[A-H]\d\d$//; push @{$family2ind{"parent2"}{$FamilyName}}, $Parent2;}
					}				
					push @{$family2ind{"progeny"}{$FamilyName}}, $IndividualName;
				}
				else
				{
					print "In excel family sheet, individual '$IndividualName' does not have family name value\n";
					exit;
				}
			}

		}

	}
	elsif ($name=~/family_trait/i)  
	{
		unless (( $data[0][0] =~/Family/i) && ( $data[0][1] =~/Trait/i)  && ( $data[0][2] =~/Donor/i) )
		{
			print "Error: In the sheet 'family_trait', the top rows must have the value 'FamilyName' 'Trait' 'Donor'\n";
			exit;
		}

		for (my $row=1; $row<=$#data; $row++)
		{
			my ($FamilyName, $Trait, $Donor) = @{$data[$row]};
			

			if ((defined $FamilyName) && ($FamilyName=~/\w/) )
			{
				$FamilyName=~s/\s//g; 
				$Trait=~s/\s//g;

				if ($Trait=~/\w/) 
				{
					unless (exists $family2traits{$FamilyName}) 
					{
						push @families, $FamilyName;
					}
					unless (defined $Donor) {$Donor= "";}
					push @{$family2traits{$FamilyName}}, "$Trait\t$Donor";
				}
				else
				{
					print "In excel family_trait sheet, family '$FamilyName' does not have trait value\n";
					exit;
				}
			}

		}

	}
	elsif ($name=~/Seq_file/i)  
	{
		unless (( $data[0][0] =~/IndividualName/i) && ( $data[0][1] =~/File1/i)  && ( $data[0][2] =~/File2/i) && ( $data[0][3] =~/Directory/i))
		{
			print "Error: In the sheet 'Seq_file', the top rows must have the value 'IndividualName' 'File1' 'File2' 'Directory' \n";
			exit;
		}

		for (my $row=1; $row<=$#data; $row++)
		{
			my ($IndividualName, $File1, $File2, $Directory) = @{$data[$row]};

			if ((defined $IndividualName) && ($IndividualName=~/\w/) )
			{
				if (($Directory=~/\w/) && ( $File1=~/\w/))
				{
					unless (-e "$Directory/$File1") 
					{
						print "$Directory/$File1 does not exist\n";
					    exit;
					}
					if ((defined $File2) && ($File2=~/\w/)) 
					{
						unless (-e "$Directory/$File2") 
						{
							print "$Directory/$File2 does not exist\n";
							exit;
						}
					}
					else
					{
						$File2= "";
					}
					$IndividualName=~s/\s//g; 
					$File1=~s/\s//g;
					$File2=~s/\s//g;
					$Directory=~s/\s//g;

					push @{$ind2file{$IndividualName}}, "$Directory\t$File1\t$File2";
				}
				else
				{
					print "In excel Seq_file sheet, individual '$IndividualName' does not have directory or file1 value\n";
					exit;
				}
			}

		}
	}

}

sub process_fastq
{
	my $fastq =shift;
	my $matchstrs = shift;
	my @matchstrs = ();
	foreach  (@{$matchstrs}) 
	{
		my @t = /(\w+)/g;
		push @matchstrs, (join '|', @t);

	}
	if ($fastq=~/gz$/) 
	{
		open ( IN,  "gunzip -c $fastq|");
	}
	else
	{
		open ( IN,  $fastq);
	}
	my $total=0;
	
	my @matched;
	for (my $i=0; $i<=$#matchstrs; $i++) 
	{
		$matched[$i]= 0;
	}

	LOOP3:while (<IN>) 
	{
		my $seq =<IN>;
		<IN>; 
		my $qual=<IN>;
		my @ascii = grep { $_ < $minQ } unpack("C*", $qual);
		if (@ascii>3) 
		{
			next LOOP3;
		}
		$total ++;

		for (my $i=0; $i<=$#matchstrs; $i++) 
		{
			my $str= $matchstrs[$i];
			if ($seq=~/$str/) 
			{
				$matched[$i]++;
			}
		}
	}
	return $total, \@matched;
}

sub process_fastq_ac
{
	my $fastqfile = shift;
	my %marker2allelecount;
	my $totalcount = 1;

	$totalcount = `wc -l $fastqfile`;
	$totalcount =~s/\D//g;

	foreach my $marker(@allmarkers)
	{
		my ($primer1, $primer2) = split "\t", $marker2primers{$marker};
		$primer2 = revcom($primer2);
		my $exec_cutadapt_cmd = $cutadapt_cmd;
		$exec_cutadapt_cmd=~s/XXX/$fastqfile/;
		$exec_cutadapt_cmd=~s/RRR/$primer2/;
		$exec_cutadapt_cmd=~s/FFF/$primer1/;

		#print $exec_cutadapt_cmd, "\n";
		system( $exec_cutadapt_cmd); 
		open FASTQ, "$tmpdir/trimmed.fastq";
		my $c=0;
		while (<FASTQ>) 
		{ 
			my $seq = <FASTQ>; <FASTQ>;<FASTQ>;
			
			chomp $seq;
			my $id = $alleleHash{$seq};
			if (defined $id) 
			{
				my ($marker, $alleleid) = split "\t", $id;
				$marker2allelecount{$marker}{$alleleid} ++;
				$c++;
			}
		}
		close FASTQ;
		print "$marker Reads total: $totalcount used: $c\n";
	}

	return ($totalcount, \%marker2allelecount);
}

sub revcom
{
        my $seq = shift;
		$seq = uc $seq;
        my $seqout = "";
        my @nt = reverse ($seq=~/(\S)/g);
        foreach my $n (@nt)
        {
			if ($compnt{$n}) 
			{
				$seqout .=$compnt{$n};
			}
             
        }
        return $seqout;
}


sub print_help
{
	print "Usage: excel_amplicon.pl -i mySampleFile -o myKeyFile -m d\n";
	print "Options:\n";
	print "-h: help\n";
	print "-i: input excel file, must be .xlsx file\n";
	print "-a: input haplotype allele file, in 3-column tab delimited text format. 1. marker;2. allele id; 3. allele sequence\n";
	print "-o: output excel file\n";	
	print "-r: output directory for tag sequences\n";
	print "-t: temporary directory. default /tmp\n";
	print "-x: max mismatch count\n";	
	print "-s: minimum fragment size in discovery mode, default 100\n"; 
	print "-m: mode. d for discovery with no known alleles; p for production with desired alleles; c for counting known alleles \n";
}

