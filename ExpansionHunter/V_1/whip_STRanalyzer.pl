open FILE, $ARGV[0];
while (my $line = <FILE>)
{
    chomp $line;
    my @Line = split /\t/, $line;
    $sample = $Line[0];
    $sex = $Line[1];
    $cmd = "python ";
    $cmd .= "/NAS/data/personal/june_dev/Pipelines/SV_ExpansionHunter/V_1/analyze_STR.py ";
    $cmd .= "--config /NAS/data/personal/june_dev/Pipelines/SV_ExpansionHunter/V_1/config.str.yaml ";
    $cmd .= "--id $sample ";
    $cmd .= "--sex $sex ";
    $cmd .= "--bamfile /NAS/data/Analysis/$sample/$sample.recal.bam ";
    $cmd .= "--genome_build hg19 ";
    $cmd .= "--outdir_vcf /NAS/data/personal/june_dev/Pipelines/SV_ExpansionHunter/V_1/Test/$sample ";
    $cmd .= "--outdir_reviewer /NAS/data/SV/Upload/ExpansionViewer/$sample ";
    system($cmd);


    $cmd = "python ";
    $cmd .= "/NAS/data/personal/june_dev/Pipelines/SV_ExpansionHunter/V_1/tucuxi_uploader_modified.py ";
    $cmd .= "--sample $sample ";
    $cmd .= "--result /NAS/data/personal/june_dev/Pipelines/SV_ExpansionHunter/V_1/Test/$sample/$sample.expansion.annotated.json ";
    $cmd .= "--config /NAS/data/personal/june_dev/Pipelines/SV_ExpansionHunter/V_1/config.tucuxi.yaml";
    print "$cmd\n";
    system($cmd);
}