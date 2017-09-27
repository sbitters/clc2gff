#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

# 2017-08-20
# STB

import sys
import regex as re
import pandas as pd
from argparse import ArgumentParser
from file_read_backwards import FileReadBackwards


def parse_input():

    parser = ArgumentParser()
    parser.add_argument("-i", type=str, dest="gff_path")
    parser.add_argument("-o", type=str, dest="output")

    parserargs = parser.parse_args()

    try:
        gff_path = parserargs.gff_path
        output_path = parserargs.output

        if gff_path is None and output_path is None:
            parser.print_help()
            raise SystemExit

        else:
            return gff_path, output_path

    except SystemExit:
        sys.exit()


def read_and_mod_gff(annot_gff):
    # Read the GFF file and replace the values "attributes" column with just the IDs of the annotated elements
    # i.e. something like "ID=id364475;Parent=gene41724;Dbxref=GeneID:19989172;exon_number=1;gbkey=exon;gene=cox2"
    # becomes "id364475"

    clc_name_corr = re.compile(r"(?<=Name=Gene) (?=[0-9]+(\n|;))")
    clc_gene_corr = re.compile(r"(?<=gene=Gene) (?=[0-9]+(\n|;))")

    id_regex = re.compile(r"(?<=ID=).+?(?=(\n|;))")
    parent_regex = re.compile(r"(?<=Parent=).+?(?=(\n|;))")
    gene_regex = re.compile(r"(?<=gene=).+?(?=(\n|;))")
    name_regex = re.compile(r"(?<=Name=).+?(?=(\n|;))")

    id_dict = {}
    gene_id_dict = {}

    id_convert = {"mRNA": "rna", "RNA": "rna", "ncRNA": "rna", "lncRNA": "rna",
                  "exon": "id", "region": "id", "CDS": "cds", "gene": "gene",
                  "variant": "var"}
    last_name = ""

    gff_line_list = []
    print("Reading GFF...")
    gff_in = FileReadBackwards(annot_gff, encoding="utf-8")
    for line in gff_in:

        if re.match("\w", line) and not line.startswith('#'):
            tab_elements = line.split("\t")

            if tab_elements[2] != "CDS":
                tab_elements[2] = tab_elements[2][0].lower() + tab_elements[2][1::]

            type = tab_elements[2]

            try:
                attributes = re.sub(clc_name_corr, "", tab_elements[-1])
            except:
                attributes = tab_elements[-1]

            try:
                attributes = re.sub(clc_gene_corr, "", attributes)
            except:
                pass

            try:
                element_id = re.search(id_regex, attributes).group()

                if element_id.startswith("mRNA"):
                    element_id = element_id[1::].lower()
            except:
                element_id = "."

            if element_id == ".":
                id_dict[id_convert[type]] = id_dict.get(id_convert[type], -1) + 1
                element_id = id_convert[type] + str(id_dict[id_convert[type]])

            try:
                gene_sym = re.search(gene_regex, attributes).group()
            except:
                gene_sym = "."

            try:
                gene_name = re.search(name_regex, attributes).group()
                if type == "gene":
                    gene_id_dict[gene_name] = element_id

                last_name = gene_name
            except:
                gene_name = "."

            if gene_name == ".":
                gene_name = last_name

            if gene_name != "." and gene_sym == ".":
                gene_sym = gene_name

            try:
                parent_gene = re.search(parent_regex, attributes).group()
            except:
                parent_gene = "."

            if parent_gene == "." and gene_name != ".":
                try:
                    parent_gene = gene_id_dict[gene_name]
                except KeyError:
                    parent_gene = "gene" + str(int(id_dict["gene"]) + 1)

            attributes = gene_name + ";" + gene_sym + ";" + parent_gene

            tab_elements[-1] = attributes

            gff_line_list.append(tab_elements)

    const = 1000000

    id_dict = {}

    last_parent = ""
    last_rna = ""

    gff_line_mod_list = []
    newgene_line_list = []
    chrom = 1

    for line in gff_line_list[::-1]:
        type = line[2]
        attributes = line[-1]

        gene_name, gene_sym, parent_gene = attributes.split(";")

        id_dict[id_convert[type]] = int(id_dict.get(id_convert[type], const)) + 1
        element_id = id_convert[type] + str(id_dict[id_convert[type]])

        if type == "mRNA":
            last_rna = element_id
            last_parent = "gene" + str(int(id_dict.get("gene", const)) + 1)
        if type == "gene" and gene_sym.startswith("Gene"):
            element_id = last_parent
        else:
            parent_gene = last_rna

        if type == "gene":
            attributes = "ID=" + element_id + ";" + "Name=" + gene_name + ";" + "gene=" + gene_sym + ";" + "gbkey=Gene"
        elif type == "mRNA":
            attributes = "ID=" + element_id + ";" + "Parent=" + "gene" + str(int(id_dict.get("gene", const)) + 1) + ";" + "gene=" + gene_sym + ";" + "gbkey=mRNA"
        elif type == "exon":
            attributes = "ID=" + element_id + ";" + "Parent=" + parent_gene + ";" + "gene=" + gene_sym
        elif type == "CDS":
            attributes = "ID=" + element_id + ";" + "Parent=" + parent_gene + ";" + "gene=" + gene_sym + ";" + "gbkey=CDS"
        elif type == "region":
            attributes = "ID=" + element_id + ";" + "Dbxref=taxon:4072;" + "chromosome=" + str(chrom)
            chrom += 1
        elif re.match("RNA", type):
            attributes = "ID=" + element_id + ";" + "Parent=" + parent_gene + ";" + "gene=" + gene_sym

        line[-1] = attributes

        gff_line_mod_list.append(line)

        if re.match("Gene", gene_sym):
            newgene_line_list.append(line)


    gff_cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    gff_df = pd.DataFrame(gff_line_mod_list, columns=gff_cols)

    newgene_df = pd.DataFrame(newgene_line_list, columns=gff_cols)

    return gff_df, newgene_df


def main():
    gff_path, output_path = parse_input()

    gff_df, newgene_df = read_and_mod_gff(gff_path)

    print("Saving GFF3...")
    gff_df.to_csv(output_path, sep="\t", index=False, header=False)
    newgene_df.to_csv((output_path[0:-5]+"_new.gff3"), sep="\t", index=False, header=False)
main()
