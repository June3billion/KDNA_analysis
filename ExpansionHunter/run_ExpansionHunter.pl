my $tool_expansionhunter = "/NAS/data/etc/ExpansionHunter/build/install/bin/ExpansionHunter";
my $tool_reviewer = "/NAS/data/etc/REViewer/build/install/bin/REViewer";
my $tool_samtools = "/NAS/data/etc/samtools-1.9/bin/samtools";
my $tool_qsub = "/NAS/sge/bin/lx-amd64/qsub";

my $file_samplelist = "/GDS/data/personal/june/Analysis/KDNA/ExpansionHunter/list_gebraID_and_sex.txt";
my $file_reference = "/NAS/data/etc/reference/hg38_mainchr/Homo_sapiens_assembly38.main_chr.fasta";
#my $file_catalog = "/NAS/data/etc/ExpansionHunter/variant_catalog/hg38/variant_catalog_updated_withofftarget_hg38.json";
my $file_catalog = "/NAS/data/etc/ExpansionHunter/variant_catalog/hg38/variant_catalog_updated_hg38.json";

my $dir_working = "/data/Analysis/Temp/";
my $dir_qsub = "/GDS/data/personal/june/Analysis/KDNA/ExpansionHunter/Qsub/";
my $dir_input = "/GDS/data/personal/june/Analysis/KDNA/Resource/GebraID/";
my $dir_output  = "/GDS/data/personal/june/Analysis/KDNA/ExpansionHunter/Output/";


my $option_thread = 3; 
my $option_locus = "AFF2,AR,ARX_1,ARX_2,ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,";
$option_locus .= "ATXN8OS,BEAN1,C9ORF72,CACNA1A,CBL,CNBP,COMP,DAB1,DIP2B,DMD,";
$option_locus .= "DMPK,FMR1,FOXL2,FXN,GIPC1,GLS,HOXA13_1,HOXA13_2,HOXA13_3,HOXD13,";
$option_locus .= "HTT,JPH3,LRP12,MARCHF6,NIPA1,NOP56,NOTCH2NLC,NUTM2B-AS1,PABPN1,PHOX2B,";
$option_locus .= "PPP2R2B,PRDM12,RAPGEF2,RFC1,RUNX2,SAMD12,SOX3,STARD7,TBP,TBX1,";
$option_locus .= "TCF4,TNRC6A,XYLT1,YEATS2,ZIC2,ZIC3";
my $option_qsub = "-S /bin/bash -l hostname=adenine21 -pe smp $option_thread"; 

open FILE, $file_samplelist;
while(my $line = <FILE>)
{
	chomp $line;
	my @Line = split /\t/, $line;
	my $sampleID = $Line[0];
	my $sex = $Line[1];
	my $file_qsub = "$dir_qsub/STR_$sampleID.sh";
	my $file_bam_input = "$dir_working/STR_$sampleID/$sampleID.hg38.recal.bam";
	my $file_bam_realigned = "$dir_working/STR_$sampleID/$sampleID.expansions_realigned.bam";
	my $file_bam_sorted = "$dir_working/STR_$sampleID/sorted_$sampleID.expansions_realigned.bam";
	my $file_vcf_str = "$dir_working/STR_$sampleID/$sample.expansions.vcf";
	my $dir_reviewer_result = "/NAS/data/SV/Upload/ExpansionViewer/$sampleID/";

	open QSUB, ">$file_qsub";
	print QSUB "date \n\n";

	print QSUB "mkdir $dir_working/STR_$sampleID/ \n";
	#print QSUB "rsync -aiL $dir_input/$sampleID* $dir_working/STR_$sampleID/ \n";
	print QSUB "date \n\n";

	print QSUB "$tool_expansionhunter \\\n";
	#print QSUB " --reads $file_bam_input \\\n";
	print QSUB " --reads $dir_input/$sampleID.hg38.recal.bam \\\n";
	print QSUB " --reference $file_reference \\\n";
	print QSUB " --variant-catalog $file_catalog \\\n";
	print QSUB " --output-prefix $dir_working/STR_$sampleID/$sampleID.expansions \\\n";
	print QSUB " --sex $sex \\\n";
	print QSUB " --threads $option_thread \n";
	print QSUB "date \n\n";
	
	print QSUB "$tool_samtools sort \\\n";
	print QSUB " $file_bam_realigned \\\n";
	print QSUB " -o $file_bam_sorted \n";
	print QSUB "$tool_samtools index \\\n";
	print QSUB " $file_bam_sorted\n";
	print QSUB "date  \n\n";

	print QSUB "rm -r $dir_reviewer_result\n";
	print QSUB "mkdir $dir_reviewer_result\n";
	print QSUB "$tool_reviewer \\\n";
	print QSUB " --reads $file_bam_sorted \\\n";
	print QSUB " --vcf $file_vcf_str \\\n";
	print QSUB " --reference $file_reference \\\n";
	print QSUB " --catalog $file_catalog \\\n";
	print QSUB " --locus $option_locus \\\n";
	print QSUB " --output-prefix $dir_reviewer_result/$sampleID \n";
	print QSUB "date\n\n";

	print QSUB "rm $file_bam_input\n";
	print QSUB "rm $file_bam_realigned\n";
	print QSUB "mkdir $dir_output/$sampleID/\n";
	print QSUB "rsync -aiL $dir_working/STR_$sampleID/  $dir_output/$sampleID/ \n";
	print QSUB "rm -r $dir_working/STR_$sampleID/ \n";
	print QSUB "date\n\n ";	

	close QSUB;

	$command_qsub = "$tool_qsub ";
	$command_qsub .= "$option_qsub ";
	$command_qsub .= "-e $dir_qsub -o $dir_qsub ";
	$command_qsub .= "$file_qsub";
	print "$command_qsub\n";
}
close FILE;
