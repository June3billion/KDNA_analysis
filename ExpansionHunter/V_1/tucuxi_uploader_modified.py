import argparse
import gzip
import json
import sys
from datetime import datetime
import pytz
import os

import boto3
import yaml


def parse_args():
    parser = argparse.ArgumentParser(description="## CNV conifer pipeline")
    sample_msg = "provide sample name (ex: ETA21-ABCD)"
    parser.add_argument("--sample", type=str, help=sample_msg)
    config_msg = "provide config.yaml file"
    parser.add_argument("--config", type=str, help=config_msg)
    result_json_msg = "provide result.json file"
    parser.add_argument("--result", type=str, help=result_json_msg)
    args = parser.parse_args()
    return args


class TucuxiCnvUploader:
    def __init__(self, sample: str, result_json_file: str, config_file: str):
        self.now = str(datetime.utcnow().isoformat() + "Z")
        self.sample = sample
        self.result_json = self.sort_json(
            self.load_json_result(result_json_file)
        )
        print(len(self.result_json), result_json_file)
        self.result_dir = os.path.dirname(os.path.abspath(result_json_file))
        self.result_txt = f"{self.result_dir}/{sample}.result.txt"
        self.result_txt_ori = f"{self.result_dir}/{sample}.result.txt.ori"
        self.intersect_json = (
            f"{self.result_dir}/{sample}.3bcnv.intersect.result.json"
        )
        self.config = self.parse_config(config_file)
        print("init_table")
        self.tucuxi_analysis_table = self.init_table()
        print(self.tucuxi_analysis_table)
        print("make_sig_gene")
        self.sig_gene_list = self.make_sig_gene()

    def init_table(self):
        dynamodb = boto3.resource(
            "dynamodb",
            region_name=self.config["RegionName"],
            aws_access_key_id=self.config["AwsAccessKeyId"],
            aws_secret_access_key=self.config["AwsSecretAccessKey"],
        )
        tucuxi_analysis_table = dynamodb.Table(
            self.config["TucuxiAnalysisTable"]
        )
        return tucuxi_analysis_table

    def parse_config(self, config_file: str) -> dict:
        config = dict()
        with open(config_file) as handle:
            config = yaml.load(handle, Loader=yaml.FullLoader)
        return config

    def make_gene_list_data(self, gene_list_file) -> dict:
        gene_set = set()
        gene_list_data_raw = list()
        gene_list_data = dict()
        cnt = 0
        with gzip.open(gene_list_file, "rt") as handle:
            header = handle.readline().strip().split("\t")
            gene_idx = header.index("gene_symbol")
            chrom_idx = header.index("chrom")
            start_idx = header.index("start")
            for line in handle:
                data = line.strip().split("\t")
                genes = data[gene_idx].split("::")
                chrom = int(
                    data[chrom_idx]
                    .replace("X", "23")
                    .replace("Y", "24")
                    .replace("MT", "25")
                )
                start = int(data[start_idx])
                for gene in genes:
                    if gene not in gene_set:
                        gene_list_data_raw.append(
                            {"gene": gene, "chrom": chrom, "start": start}
                        )
                        gene_set.add(gene)
        for elem in sorted(
            gene_list_data_raw, key=lambda x: (x["chrom"], x["start"])
        ):
            elem["order"] = cnt
            gene_list_data[elem["gene"]] = elem
            cnt += 1

        return gene_list_data

    def count_grey(self):
        """
        remain_cnt, filtered_cnt 로 return 한다.
        """
        with open(self.result_txt_ori) as handle:
            cnt_total = len(handle.readlines())
        with open(self.result_txt) as handle:
            cnt_filtered = len(handle.readlines())
        return cnt_filtered, cnt_total - cnt_filtered

    def sort_json(self, result_json):
        sorted_result_json = sorted(
            result_json,
            key=lambda x: (
                (
                    x["pos"]
                    .split(":")[0]
                    .replace("X", "23")
                    .replace("Y", "24")
                    .replace("M", "25")
                ),
                (x["pos"].split(":")[1].split("-")[0]),
            ),
        )
        return sorted_result_json

    def load_json_result(self, result_json_file: str) -> dict:
        result_json = dict
        with open(result_json_file) as handle:
            result_json = json.load(handle)
        return result_json

    def make_sig_gene(self):
        sig_gene_list = list()
        # sig_gene_data_raw = self.load_json_result(self.intersect_json)
        # for _, genes in sig_gene_data_raw.items():
        #   for gene in genes:
        #        sig_gene_list.append(gene)
        return sig_gene_list

    def make_gene_disease_data(self, disease_file: str):
        gene_disease_data = dict()
        with gzip.open(disease_file, "rt") as handle:
            header = handle.readline().strip().split("\t")
            omimPhenoId_idx = header.index("#omimPhenoId")
            orphaId_idx = header.index("orphaId")
            geneSymbol_idx = header.index("geneSymbol")
            inheritance_idx = header.index("inheritances:value")
            title_idx = header.index("title")
            onset_age_idx = header.index("onsetAges:value")
            hgnc_idx = header.index("hgncId")
            ncbi_idx = header.index("ncbiGeneId")
            ensembl_idx = header.index("ensemblGeneId")
            omim_gene_idx = header.index("omimGeneId")
            orpha_gene_idx = header.index("orphaGeneSymbol")

            for line in handle:
                data = line.strip().split("\t")
                omim_pheno_id = data[omimPhenoId_idx]
                orpha_id = data[orphaId_idx]
                if omim_pheno_id == "-" and orpha_id == "-":
                    continue
                gene_symbol = data[geneSymbol_idx]
                inheritance = data[inheritance_idx].split("||")
                title = data[title_idx]
                onset_age = data[onset_age_idx]
                hgnc = data[hgnc_idx]
                ncbi = data[ncbi_idx]
                ensembl = data[ensembl_idx]
                omim_gene = data[omim_gene_idx]
                orpha_gene = data[orpha_gene_idx]
                if gene_symbol not in gene_disease_data:
                    gene_disease_data[gene_symbol] = dict()
                gene_disease_data[gene_symbol][
                    f"{gene_symbol}_{omim_pheno_id}"
                ] = {
                    "inheritance": inheritance,
                    "title": title,
                    "phenoId": {"ORPHA": orpha_id, "OMIM": omim_pheno_id},
                    "onsetAge": onset_age,
                    "ncbi": ncbi,
                    "symbol": gene_symbol,
                    "orpha": orpha_gene,
                    "ensembl": ensembl,
                    "hgnc": hgnc,
                    "omim": omim_gene,
                }
        return gene_disease_data

    def get_disease(self, gene: str, gene_disease_data: dict):
        diseases = list()
        if gene not in gene_disease_data:
            return diseases
        for gene_symbol_omim_phemo_id, pheno_data in gene_disease_data[
            gene
        ].items():
            if pheno_data["phenoId"]["OMIM"] != "-":
                phenoid = f"OMIM:{pheno_data['phenoId']['OMIM']}"
            else:
                phenoid = pheno_data["phenoId"]["ORPHA"]
            diseases.append(
                {
                    "phenoId": phenoid,
                    "ip": pheno_data["inheritance"],
                    "geneId": pheno_data["omim"],
                    "title": pheno_data["title"],
                }
            )
        return diseases

    def make_cnv_entity(self, region_set: set):
        # Write CNV version
        cnv_version_data = dict()
        cnv_version_data["PK"] = self.sample
        cnv_version_data["SK"] = f"VERSION#{self.config['Version']}"
        cnv_version_data["type"] = "cnvVersion"
        cnv_version_data["createdAt"] = self.now
        cnv_version_data["version"] = f"{self.config['Version']}"
        ####
        # (
        #    cnv_version_data["remain_cnt"],
        #    cnv_version_data["filtered_cnt"],
        # ) = self.count_grey()
        cnv_version_data["remain_cnt"] = 25
        cnv_version_data["filtered_cnt"] = 0
        ####
        print(cnv_version_data)
        self.tucuxi_analysis_table.put_item(Item=cnv_version_data)

        # Write CNV elements
        # explanation = list()
        for region_elem in self.result_json:

            ####
            # print(region_elem["pos"], region_elem["pos"] == "8:7809329-8176818")
            if region_elem["pos"] not in region_set:
                continue
            cnv_data = dict()
            cnv_data["PK"] = f"{self.sample}|||CNV#{self.config['Version']}"
            cnv_data["SK"] = f"REGION#{region_elem['pos']}"
            cnv_data["cytoband"] = f"{region_elem['cyto']}"
            # for evidence, evidence_detail in region_elem["detail"][
            #    "evidences"
            # ].items():
            #    rule_data = dict()
            #    rule_data["category"] = evidence  # 1A
            #    rule_data["score"] = evidence_detail["res"]
            #    cnv_data["rules"].append(rule_data)
            rule_data = dict()
            rule_data["category"] = f"{region_elem['rules']}"
            cnv_data["rules"] = list()
            cnv_data["rules"].append(rule_data)
            cnv_data["cnvVersion"] = self.config["Version"]
            cnv_data["createdAt"] = self.now
            cnv_data["score"] = f"{region_elem['score']}"
            cnv_data["consequence"] = f"{region_elem['status']}"  # del
            cnv_data["ACMG"] = f"{region_elem['acmg']}"  # Pathogenic
            cnv_data["pos"] = f"{region_elem['pos']}"
            cnv_data["id"] = self.sample
            cnv_data["tool"] = f"{self.config['Tool']}"
            cnv_data["chr"] = f"{region_elem['pos'].split(':')[0]}"
            # [cnv_data["start"], cnv_data["end"],] = (
            [cnv_data["start"],] = (
                region_elem["pos"].split(":")[1].split("-")
            )
            cnv_data["type"] = "CNV"
            ####
            print(cnv_data)
            ####
            self.tucuxi_analysis_table.put_item(Item=cnv_data)

    def make_cnv_gene_entity(
        self, sig_gene_list: list, gene_disease_data: dict, gene_list_data: dict
    ):
        region_set = set()
        for region_elem in self.result_json:
            # cnv_gene_data = dict()
            region_genes = region_elem["genes"]

            for gene in region_genes:
                print(gene)
                ####
                # if gene not in sig_gene_list:
                #    continue
                region_set.add(region_elem["pos"])
                upload_data = dict()
                upload_data[
                    "PK"
                ] = f"{self.sample}|||CNV#{self.config['Version']}"
                upload_data["SK"] = f"REGION#{region_elem['pos']}|||GENE#{gene}"
                upload_data["createdAt"] = self.now
                upload_data["type"] = "CNVGene"
                upload_data["id"] = self.sample
                upload_data["cnvVersion"] = self.config["Version"]
                upload_data["tool"] = "EH"
                upload_data["symbol"] = gene
                upload_data["diseases"] = self.get_disease(
                    gene, gene_disease_data
                )
                try:
                    order = gene_list_data[gene]["order"]
                except KeyError:
                    order = 0
                upload_data["order"] = order
                ####
                # print(upload_data)
                ####
                self.tucuxi_analysis_table.put_item(Item=upload_data)
        return region_set

    def main(self):
        disease_file = "/NAS/data/CNV/conifer/data/disease.txt.gz"
        gene_list_file = "/data/DB/processedData/RefSeq/result/merged.refseq.grch37.all.txt.gz"
        gene_list_data = self.make_gene_list_data(gene_list_file)
        gene_disease_data = self.make_gene_disease_data(disease_file)
        sig_gene_list = self.make_sig_gene()
        region_set = self.make_cnv_gene_entity(
            sig_gene_list, gene_disease_data, gene_list_data
        )
        # print(region_set)
        self.make_cnv_entity(region_set)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            f"#usage: python {sys.argv[0]} --sample [sample] --result [result.json] --config [config.yaml]"
        )
        print(
            f"ex) python {sys.argv[0]} --sample ETC21-PPJZ --result /data/ken_dev/pipeline/conifer_pipeline/result/210419/ETC21-PPJZ.result.acmg.json --config config.tucuxi.yaml"
        )
        print(
            f"ex) python {sys.argv[0]} --sample WILL --result /NAS/data/Analysis/WILL/CNV/WILL.result.acmg.json --config config.tucuxi.yaml"
        )

        sys.exit()

    ARGS = parse_args()
    sample = ARGS.sample
    result_json_file = ARGS.result
    config_file = ARGS.config

    tucuxi_cnv_uploader = TucuxiCnvUploader(
        sample, result_json_file, config_file
    )
    tucuxi_cnv_uploader.main()
