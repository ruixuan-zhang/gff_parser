import pandas as pd
from Bio import SeqIO
from . import parse_gff3
from collections import defaultdict
import os 

def make_cds_dict(genome_path, gff_path):
    
    # Read genome fasta into dictionary, each contig_id is a key
    handle_genome_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

    # Read genome gff into panda dataframe
    handle_gff = parse_gff3(gff_path, verbose = False, parse_attributes = True)
    handle_CDS = handle_gff.copy()
    handle_CDS = handle_CDS[handle_CDS["gbkey"] == "CDS"]
    handle_CDS = handle_CDS.reset_index(drop=True)

    CDS_strand_dict = defaultdict(str)
    CDS_seqid_dict = defaultdict(str)
    CDS_region_dict = defaultdict(list)
    CDS_misc_dict = defaultdict(list)
    
    for i in handle_CDS.index:
        # keys 
        seqid = handle_CDS.iloc[i, 0]
        strand = handle_CDS.iloc[i, 6]
        start = handle_CDS.iloc[i, 3]
        end = handle_CDS.iloc[i, 4]
        protein_id = handle_CDS.iloc[i, 9]
        gene = handle_CDS.iloc[i, 15]
        locus_tag = handle_CDS.iloc[i, 19]

        CDS_strand_dict[locus_tag] = strand
        CDS_seqid_dict[locus_tag] = seqid
        CDS_region_dict[locus_tag].append([start, end])
        CDS_misc_dict[locus_tag] = [protein_id, gene]

    return(CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, CDS_misc_dict, handle_genome_dict)

def extract_cds(CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, handle_genome_dict):
    CDS_seq_dict = defaultdict(str)
    for gene_record in CDS_seqid_dict.keys():
        seqid = CDS_seqid_dict[gene_record]
        contig_seq = handle_genome_dict[seqid].seq
        strand = CDS_strand_dict[gene_record]

        CDS_seq_list = []

        if strand == "+":    
            for region_list in CDS_region_dict[gene_record]:
                start = int(region_list[0])
                end = int(region_list[1])
                CDS_seq_list.append(str(contig_seq[start-1: end]))
        else:
            for region_list in CDS_region_dict[gene_record]:
                start = int(region_list[0])
                end = int(region_list[1])
                CDS_seq_list.append(str(contig_seq[start-1: end].reverse_complement()))
        
        CDS_seq = "".join(CDS_seq_list)
        CDS_seq_dict[gene_record] = CDS_seq
    
    return CDS_seq_dict

def extract_upstream(CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, left_shift, handle_genome_dict):
    UP_seq_dict = defaultdict(str)
    for gene_record in CDS_seqid_dict.keys():
        seqid = CDS_seqid_dict[gene_record]
        contig_seq = handle_genome_dict[seqid].seq
        strand = CDS_strand_dict[gene_record]

        if strand == "+":
            start = int(CDS_region_dict[gene_record][0][0]) - 1 
            left_bound = start + left_shift
            up_seq = str(contig_seq[left_bound: start])
        else: 
            start = int(CDS_region_dict[gene_record][-1][1])
            right_bound = start - left_shift
            up_seq = str(contig_seq[start : right_bound].reverse_complement())
        
        UP_seq_dict[gene_record] = up_seq
    
    return UP_seq_dict 

def extract_after_start(CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, right_shift, handle_genome_dict):
    DOWN_seq_dict = defaultdict(str)
    for gene_record in CDS_seqid_dict.keys():
        seqid = CDS_seqid_dict[gene_record]
        contig_seq = handle_genome_dict[seqid].seq
        strand = CDS_strand_dict[gene_record]

        if strand == "+":
            start = int(CDS_region_dict[gene_record][0][0]) - 1 # to not include start codon 
            right_bound = start + right_shift
            DOWN_seq = str(contig_seq[start: right_bound])
        else: 
            start = int(CDS_region_dict[gene_record][-1][1])
            left_bound = start - right_shift
            DOWN_seq = str(contig_seq[left_bound : start].reverse_complement())
        
        DOWN_seq_dict[gene_record] = DOWN_seq
    
    return DOWN_seq_dict 

def extract_around_start(UP_seq_dict, DOWN_seq_dict):
    AROUND_start_dict = defaultdict(str)
    for gene_record in CDS_strand_dict.keys():
        strand = CDS_strand_dict[gene_record]
        around_start_seq = UP_seq_dict[gene_record] + DOWN_seq_dict[gene_record]
        AROUND_start_dict[gene_record] = around_start_seq
        
    return AROUND_start_dict


def run_all(genome_path, gff_path, left_shift, right_shift, working_dir, seq_name):
    CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, CDS_misc_dict, handle_genome_dict = make_cds_dict(genome_path, gff_path)

    # Extract the CDS sequences
    CDS_seq_dict = extract_cds(CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, handle_genome_dict)

    # Extract the upstream sequences with a shift of 5 (for example)
    UP_seq_dict = extract_upstream(CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, left_shift, handle_genome_dict)

    # Extract the downstream sequences with a shift of 5 (for example)
    DOWN_seq_dict = extract_after_start(CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, right_shift, handle_genome_dict)
    
    AROUND_start_dict = extract_around_start(UP_seq_dict, DOWN_seq_dict)
    
    # output 
    cds_folder = os.path.join(working_dir, "cds")
    n_term_long_folder = os.path.join(working_dir, "n_term_long")
    after_start_long_folder = os.path.join(working_dir, "after_start_long")

    os.makedirs(cds_folder, exist_ok=True)
    os.makedirs(n_term_long_folder, exist_ok=True)
    os.makedirs(after_start_long_folder, exist_ok=True)

    # write fna files 
    with open(os.path.join(cds_folder, seq_name + ".fna"), "w") as cds_out:
        for key, values in CDS_seq_dict.items():
            cds_out.write(">" + key + "\n")
            cds_out.write(values + "\n")

    with open(os.path.join(n_term_long_folder, seq_name + ".fna"), "w") as n_term_out:
        for key, values in DOWN_seq_dict.items():
            n_term_out.write(">" + key + "\n")
            n_term_out.write(values + "\n")
        
    with open(os.path.join(after_start_long_folder, seq_name + ".fna"), "w") as after_out:
        for key, values in AROUND_start_dict.items():
            after_out.write(">" + key + "\n")
            after_out.write(values + "\n")