import os
import sys
import json
import yaml
import argparse
import subprocess

from readline import replace_history_item


def get_arguments() -> argparse.Namespace:
    """
    Guide and check arguments required for this script.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        help="input config file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--id",
        help="input sample ID",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--sex",
        help="input sample sex",
        choices=["male", "female"],
        type=str,
        required=True,
    )
    parser.add_argument(
        "--bamfile",
        help="input bam file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--genome_build",
        help="reference version, GRCh37 or GRCh38",
        choices=["hg19", "hg38"],
        type=str,
        required=True,
    )
    parser.add_argument(
        "--outdir_vcf",
        help="output directory for vcf and bam",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--outdir_reviewer",
        help="output directory for reviewer",
        type=str,
        required=True,
    )
    return parser.parse_args()


def run_expansionhunter(config: dict, args: argparse.Namespace):
    """Run expansion hunter for the input sample

    Note:
        Required: (Prepared with config file)
            Expansion Hunter:

    ARGS:

    Example
    """
    expansionhunter = config["TOOL"]["EXPANSIONHUNTER"]
    bamfile = args.bamfile
    sex = args.sex
    output_prefix = f"{args.outdir_vcf}/{args.id}.expansion"
    thread = config["OPTION"]["THREAD"]
    catalog = config["FILE"]["CATALOG"]["HG19"]
    reference = config["FILE"]["REFERENCE"]["HG19"]
    if args.genome_build == "hg38":
        catalog = config["FILE"]["CATALOG"]["HG38"]
        reference = config["FILE"]["REFERENCE"]["HG38"]

    cmd_expansionhunter = (
        f"{expansionhunter} "
        f"--reads {bamfile} "
        f"--sex {sex} "
        f"--reference {reference} "
        f"--variant-catalog {catalog} "
        f"--output-prefix {output_prefix} "
        f"--thread {thread} "
    )
    print(cmd_expansionhunter)
    subprocess.run(cmd_expansionhunter, shell=True, check=True)
    return


def run_reviewer(config: dict, args: argparse.Namespace):
    """Run reviewr, by locus
    ddd

    """
    reviewer = config["TOOL"]["REVIEWER"]
    samtools = config["TOOL"]["SAMTOOLS"]
    bamfile_realigned = f"{args.outdir_vcf}/{args.id}.expansion_realigned.bam"
    bamfile_sorted = f"{args.outdir_vcf}/{args.id}.expansion_sorted.bam"
    vcf_str = f"{args.outdir_vcf}/{args.id}.expansion.vcf"
    output_prefix = f"{args.outdir_reviewer}/{args.id}"
    catalog = config["FILE"]["CATALOG"]["HG19"]
    reference = config["FILE"]["REFERENCE"]["HG19"]
    if args.genome_build == "hg38":
        catalog = config["FILE"]["CATALOG"]["HG38"]
        reference = config["FILE"]["REFERENCE"]["HG38"]
    loci = config["OPTION"]["LOCI"].split(",")

    cmd_sorting = f"{samtools} sort {bamfile_realigned} -o {bamfile_sorted}"
    cmd_indexing = f"{samtools} index {bamfile_sorted}"
    print(cmd_sorting)
    print(cmd_indexing)
    subprocess.run(cmd_sorting, shell=True, check=True)
    subprocess.run(cmd_indexing, shell=True, check=True)

    for locus in loci:
        cmd_reviewer = (
            f"{reviewer} "
            f"--reads {bamfile_sorted} "
            f"--vcf {vcf_str} "
            f"--reference {reference} "
            f"--catalog {catalog} "
            f"--locus {locus} "
            f"--output-prefix {output_prefix} "
        )
        print(cmd_reviewer)
        subprocess.run(cmd_reviewer, shell=True)
    return


def annotate_str(config: dict, args: argparse.Namespace):
    """This is qwerty
    ccddd
    """
    path_vcf = f"{args.outdir_vcf}/{args.id}.expansion.vcf"
    path_annotated = f"{args.outdir_vcf}/{args.id}.expansion.annotated.json"
    path_annotated_temp = (
        f"{args.outdir_vcf}/{args.id}.3bcnv.intersect.result.json"
    )
    str_db = config["FILE"]["ANNOTATION"]
    file_vcf = open(path_vcf, "r")
    vcf_lines = file_vcf.readlines()
    file_vcf.close()
    file_str_db = open(str_db, "r")
    lines_str_db = file_str_db.readlines()
    file_str_db.close()

    output = []

    for vcf_line in vcf_lines:
        vcf_line = vcf_line.rstrip("\n")
        if vcf_line.startswith("#"):
            continue
        # print(vcf_line)
        (
            chromosome,
            position,
            id,
            ref,
            alt,
            quality,
            filter,
            info,
            format,
            calls,
        ) = vcf_line.split("\t")
        chromosome_noprefix = chromosome.replace("chr", "")
        (
            end,
            end_val,
            ref,
            ref_val,
            rl,
            rl_val,
            ru,
            ru_val,
            varid,
            varid_val,
            repid,
            repid_val,
        ) = info.replace("=", ";").split(";")

        (
            gt,
            so,
            repcn,
            repci,
            adsp,
            adfl,
            adir,
            lc,
        ) = calls.split(":")
        max_copy = max(repcn.split("/"))
        if filter != "PASS":
            continue

        for line_str_db in lines_str_db:
            line_str_db = line_str_db.rstrip("\n")
            (
                locus,
                gene_omim,
                disorder,
                disorder_omim,
                onset,
                repeat,
                inheritance,
                normal,
                pathogenic,
            ) = line_str_db.split("\t")
            max_normal = max(normal.split("-"))
            min_pathogenic = min(pathogenic.split("-"))
            gene = locus.split("_")[0]
            if locus != varid_val:
                continue
            if int(max_copy) >= int(min_pathogenic):
                output.append(
                    dict(
                        acmg="Pathogenic",
                        cyto=f"http://bamvi.3billion.io/expansion_review?sample={args.id}&gene={locus},__  __,https://stripy.org/database/{locus}",
                        detail="NA",
                        genes=[gene],
                        num_genes=1,
                        pos=f"{chromosome_noprefix}:{position},{ru_val}",
                        rules=f"{repcn}({repci})",
                        score=f"{pathogenic}({normal})",
                        status="EXP",
                        sym_res=[],
                    )
                )
                continue

            if int(max_copy) > int(max_normal):
                output.append(
                    dict(
                        acmg="VUS",
                        cyto=f"http://bamvi.3billion.io/expansion_review?sample={args.id}&gene={locus},__  __,https://stripy.org/database/{locus}",
                        detail="NA",
                        genes=[gene],
                        num_genes=1,
                        pos=f"{chromosome_noprefix}:{position},{ru_val}",
                        rules=f"{repcn}({repci})",
                        score=f"{pathogenic}({normal})",
                        status="EXP",
                        sym_res=[],
                    )
                )
                continue

    with open(path_annotated, "w") as fh:
        json.dump(output, fh, indent=4)

    with open(path_annotated_temp, "w") as fh:
        json.dump(output, fh, indent=4)


def main():
    """
    Note:
        Required: (Prepared with config file)
            ExpansionHunter
            REViwer
    Author: June (Youngjune Bhak)
    Contact: youngjune29bhak@gmail.com
    Date (ver): 2021.09.29
    """
    args = get_arguments()
    with open(args.config, "r") as fh:
        config = yaml.load(fh, Loader=yaml.FullLoader)
    if not os.path.exists(args.outdir_vcf):
        cmd_mkdir = f"mkdir -p {args.outdir_vcf}"
        print(cmd_mkdir)
        subprocess.run(cmd_mkdir, shell=True, check=True)
    if not os.path.exists(args.outdir_reviewer):
        cmd_mkdir = f"mkdir -p {args.outdir_reviewer}"
        print(cmd_mkdir)
        subprocess.run(cmd_mkdir, shell=True, check=True)
    run_expansionhunter(config, args)
    run_reviewer(config, args)
    annotate_str(config, args)


if __name__ == "__main__":
    main()
