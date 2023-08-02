import pandas as pd
from Bio import SeqIO
from . import parse_gff3
from collections import defaultdict
import os 

class CDS:
    """ 
    This class saves the coding sequence (CDS) features of a gene including 
    - seqid, the contig id in the genome fasta 
    - strand
    - region, a list of cds regions such as [[exon1_start, exon1_end] , [exon2_start, exon2_end]]
    - misc: some other features scu has the protein_id, and gene id 
    """
    def __init__(self, seqid, strand, region, misc):
        self.seqid = seqid
        self.strand = strand
        self.region = region
        self.misc = misc
        self.length = sum([end - start + 1 for start, end in region])

def make_cds_object(genome_path, gff_path):
    """ 
    This function reads genome fasts into dictionary using Bio.SeqIO, each contig_id is a key
    This function reads genome gff into panda dataframe 

    Args:
        genome_path : A string that specifies the path to the genome file
        gff_path : A string that specifies the path to the GFF file
    
    Returns:
        A CDS class object
        A dictionary of CDS
    """

    # Read genome fasta into dictionary, each contig_id is a key
    handle_genome_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

    # Read genome gff into panda dataframe
    handle_gff = parse_gff3(gff_path, verbose = False, parse_attributes = True)
    handle_CDS = handle_gff.copy()
    handle_CDS = handle_CDS[handle_CDS["gbkey"] == "CDS"]
    handle_CDS = handle_CDS.reset_index(drop=True)

    cds_objects_dict = {}

    # CDS_strand_dict = defaultdict(str)
    # CDS_seqid_dict = defaultdict(str)
    # CDS_region_dict = defaultdict(list)
    # CDS_misc_dict = defaultdict(list)
    
    for i in handle_CDS.index:
        # keys 
        seqid = handle_CDS.loc[i, "Seqid"]
        strand = handle_CDS.loc[i, "Strand"]
        start = int(handle_CDS.loc[i, "Start"])
        end = int(handle_CDS.loc[i, "End"])
        protein_id = handle_CDS.loc[i, "protein_id"]
        gene = handle_CDS.loc[i, "gene"]
        locus_tag = handle_CDS.loc[i, "locus_tag"]

        region = [[start, end]]
        misc = [protein_id, gene]
        
        if locus_tag in cds_objects_dict: 
            cds_objects_dict[locus_tag].region.append(region[0])
            cds_objects_dict[locus_tag].length += region[0][1] - region[0][0] + 1  # Update length
        else:
            cds = CDS(seqid, strand, region, misc)
            cds_objects_dict[locus_tag] = cds

        # CDS_strand_dict[locus_tag] = strand
        # CDS_seqid_dict[locus_tag] = seqid
        # CDS_region_dict[locus_tag].append([start, end])
        # CDS_misc_dict[locus_tag] = [protein_id, gene]
    
    return cds_objects_dict, handle_genome_dict
    # return(CDS_strand_dict, CDS_seqid_dict, CDS_region_dict, CDS_misc_dict, handle_genome_dict)

def extract_cds(cds_objects_dict, handle_genome_dict):
    """
    This function takes the cds_objects_dict from previous function {locus_tag : cds object}
    """

    CDS_seq_dict = defaultdict(str)
    for gene_record, cds in cds_objects_dict.items():
        contig_seq  = handle_genome_dict[cds.seqid].seq

        CDS_seq_list = []

        if cds.strand == "+":    
            for region_list in cds.region:
                start = int(region_list[0])
                end = int(region_list[1])
                CDS_seq_list.append(str(contig_seq[start-1: end]))
        else:
            for region_list in cds.region:
                start = int(region_list[0])
                end = int(region_list[1])
                CDS_seq_list.append(str(contig_seq[start-1: end].reverse_complement()))
        
        CDS_seq = "".join(CDS_seq_list)
        CDS_seq_dict[gene_record] = CDS_seq
    
    return CDS_seq_dict

def extract_upstream(cds_objects_dict, left_shift, handle_genome_dict):
    """
    This function takes in the cds_objects_dict as well and use handle_genome_dict to get the contig information

    Args: 
        cds_objects_dict
        left_shift: the upstream region length you want to extract. A negative number !!!! such as -4, -60, etc 
        handle_genome_dict

    Return
        UP_seq_dict: a dictionary saves {locus_tag : sequences }
        UP_error_dict : a dictioanry saves {locus_tag : sequences } # a common problem is that the upstream is too short to be extracted, these sequences will be saved separatedly 

    Note: 
        Due to the python list index, if you want to extract region 3 nt prior to the start codon (0), the region should be [-3, 0]
    """

    UP_seq_dict = defaultdict(str)
    UP_error_dict = defaultdict(str)

    for gene_record, cds in cds_objects_dict.items():
        contig_seq  = handle_genome_dict[cds.seqid].seq
        up_seq = ""
        up_error_seq = ""

        if cds.strand == "+":
            start = int(cds.region[0][0]) - 1 
            if start + left_shift > 0:
                left_bound = start + left_shift
                up_seq = str(contig_seq[left_bound: start])
                UP_seq_dict[gene_record] = up_seq
            else:
                up_error_seq = str(contig_seq[0: start])
                UP_error_dict[gene_record] = up_error_seq
        else: 
            start = int(cds.region[-1][1])
            right_bound = start - left_shift
            if right_bound <= len(contig_seq):
                up_seq = str(contig_seq[start : right_bound].reverse_complement())
                UP_seq_dict[gene_record] = up_seq
            else:
                up_error_seq = str(contig_seq[start: right_bound].reverse_complement())
                UP_error_dict[gene_record] = up_error_seq

        # UP_seq_dict[gene_record] = up_seq
        # UP_error_dict[gene_record] = up_error_seq
    
    return UP_seq_dict, UP_error_dict

def extract_after_start(cds_objects_dict, CDS_seq_dict, right_shift, include_start=True):
    """
    This function aims to extract the region after the start codon (including start codon )
    
    Args:
        cds_objects_dict
        right_shfit
        handle_genome_dict
        include_start

    Return: 
        DOWN_seq_dict : downstream sequences after start codon 
        DOWN_error_dict : downstream sequences if the right_shift longer than the cds length, to avoid of including stop codon, only region before stop codon will be output

    """
    DOWN_seq_dict = defaultdict(str)
    DOWN_error_dict = defaultdict(str)

    for gene_record, cds in cds_objects_dict.items():
        CDS_seq = CDS_seq_dict[gene_record]
        down_seq = ""
        down_error_seq = ""

        if include_start == True:
            if right_shift > len(CDS_seq) - 3:
                down_error_seq = str(CDS_seq[0:-3])
                DOWN_error_dict[gene_record] = down_error_seq
            else:
                down_seq = str(CDS_seq[0: right_shift])
                DOWN_seq_dict[gene_record] = down_seq
        else:
            if right_shift > len(CDS_seq) - 6: # extracted region at most will be ATG (start) | -> AAA ...... GGG <- | TGA (stop)
                down_error_seq = str(CDS_seq[3: -3])
                DOWN_error_dict[gene_record] = down_error_seq
            else:
                down_seq = str(CDS_seq[3: 3 + right_shift])
                DOWN_seq_dict[gene_record] = down_seq


    return DOWN_seq_dict, DOWN_error_dict
 

def check_region_status(UP_seq_dict, DOWN_seq_dict):
    """
    This functions aims to get which gene meet the criteria of region selection
    """

    UP_set = set(list(UP_seq_dict.keys()))
    DOWN_set = set(list(DOWN_seq_dict.keys()))
    BOTH_set = UP_set.intersection(DOWN_set)

    return BOTH_set
        
def extract_around_start(UP_seq_dict, DOWN_seq_dict, BOTH_set):
    """
    This functions aims to merge region before CDS and after CDS together
    """

    AROUND_start_dict = defaultdict(str)

    for gene_record in BOTH_set:        
        around_start_seq = UP_seq_dict[gene_record] + DOWN_seq_dict[gene_record]
        AROUND_start_dict[gene_record] = around_start_seq
        
    return AROUND_start_dict

