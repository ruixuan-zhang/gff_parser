from . import cds_extract
 
def around_start(genome_path, gff_path, left_shift, right_shift, include_start=True):
    # create CDS object
    cds_objects_dict, handle_genome_dict = cds_extract.make_cds_object(genome_path, gff_path)

    CDS_seq_dict = cds_extract.extract_cds(cds_objects_dict, handle_genome_dict)

    UP_seq_dict, UP_error_dict = cds_extract.extract_upstream(cds_objects_dict, left_shift, handle_genome_dict)

    DOWN_seq_dict, DOWN_error_dict = cds_extract.extract_after_start(cds_objects_dict, CDS_seq_dict, right_shift, include_start)

    # check which genes meet the criteria for region selection
    BOTH_set = cds_extract.check_region_status(UP_seq_dict, DOWN_seq_dict)
    # combine the upstream and downstream regions
    AROUND_start_dict = cds_extract.extract_around_start(UP_seq_dict, DOWN_seq_dict, BOTH_set)

    return AROUND_start_dict, UP_error_dict, DOWN_error_dict
