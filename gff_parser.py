import pandas as pd
import numpy as np
from tqdm import tqdm


# 处理 Attribute column 的问题: 不同 type 的 record 的 attribute 是不一样的
# function 1: Get first 8 columns (usual columns) and get attribute column
# function 2: Turn attribute column into a dictionary 
# function 3: 如果存在，就加入，没有就 Nan

def _parse_cols(lines, verbose, return_attrib_dict=False):
    """
    Hidden function to parse the tabular part of the data
    :param lines: (list of strings) - the lines of the file
    :return: (pd.DataFrame) - columns = Seqid, Source, Type, Start, End, Score, Strand, Phase

    """

    # define the name of the first 8 columns

    data = {"Seqid": [],
            "Source": [],
            "Type": [],
            "Start": [],
            "End": [],
            "Score": [],
            "Strand": [],
            "Phase": []}

    # save the attribute content in a list named as `attrib` 
    attribs_list = []
    for line in (tqdm(lines) if verbose else lines):
        if line[0] != '#':
            columns = line.split('\t')

            for i, key in enumerate(data):
                data[key].append(columns[i])

            attribs_list.append(columns[-1])


    df = pd.DataFrame(data).replace('.', np.nan)

    if return_attrib_dict:
        return df, attribs_list

    else:
        return df

def _make_dict_from_attrib_list(attribs_element):
    """
    Take the string type variable from the `attribs_list` generated in _parse_cols
    Convert the semi-colon separated key-value pairs into a dictionary for each row in the attribute column.
    """
    keys = [x.split("=")[0] for x in attribs_element.split(";")]
    vals = [x.split("=")[1] for x in attribs_element.split(";")]

    # save the keys, vals pairs in a dictionary named data 
    attrib_pair_dict = {}
    for k, v in zip(keys, vals):
        attrib_pair_dict[k] = v

    return attrib_pair_dict


def _get_attrib_keys(attrib_col, verbose):
    """
    Aim: get the all, unique keys from the input attribute list. 
    :param attrib_col: the returned list from _parse_cols function 
    :param verbose: verbose, default true
	"""

    if verbose:
        print("Finding unique attribute keys...")
    
    # define a set, keys in the attribute list will be saved here, such as ID, locus_tag, ... etc
    attrib_keys = set()

    for line in (tqdm(attrib_col) if verbose else attrib_col):
        key_val_arr = line.split(";")
        for key_val in key_val_arr:
            key, val = key_val.split("=")
            attrib_keys.add(key)

    return list(attrib_keys)

def _make_attrib_table(attribs_list, attrib_keys, verbose):
    """
    Constructs a DataFrame from the attributes column by splitting the key-value pairs.
	:param attribs_list: the returned list from _parse_cols
    :param attrib_keys: the set of attribute keys returned from _get_attrib_keys
    """
    if verbose:
        print("Making attribute table...")

    data = {}
    for key in attrib_keys:
        data[key] = []

    for string in (tqdm(attribs_list) if verbose else attribs_list):
        row_attrib_dict = _make_dict_from_attrib_list(string)

        for key in data:
            if key in row_attrib_dict:
                data[key].append(row_attrib_dict[key])

            else:
                data[key].append(np.nan)

    df = pd.DataFrame(data)
    return df

def _parse_all(lines, verbose):
    if verbose:
        print('Building structured data...')
    
    df, attribs = _parse_cols(lines, verbose, return_attrib_dict=True)
    # retrun attribs, a column that each element is a attribute record

    if verbose:
        print("Adding Supplemental Attribute table...")
    
    attrib_keys = _get_attrib_keys(attribs, verbose)
    # _get_attrib_keys return [keys] (ID, dbxref, ...)
    supplement_table = _make_attrib_table(attribs, attrib_keys, verbose)
    new_df = pd.merge(df, supplement_table, how='left', left_index=True, right_index=True)
    return new_df

def parse_gff3(filepath, verbose=True, parse_attributes=False):
    """
    This is the only public function of this package and it just parses a gff3 file into a pandas dataframe.

    **Background**  - There are 8 columns (Seqid, Source, Type, Start, End, Score, Strand, Phase) that are in every gff3 file (required by the spec).
    The last column stores additional data called 'attributes' which are stored in the last column in a weird form that is similar to key value pair tuples.
    The logic of this function is that if parameter parse_attributes is true, this will collect all this data and add it to the dataframe by creating a column
    for every unique key in all the rows. Many of these keys are always included (of at lease often) but others not so much so you can end up with some
    sparse columns occasionally.

    :param filepath: (pathlike, str) - the path to a gff3 file
    :param verbose: (bool, default = True) - This function often takes a while to run so this function prints updates and progress bars to the screen
    :param parse_attributes: (bool, default = False) - Read background in the docstring. Provides additional information but increased runtime.
    :return: pd.DataFrame
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    lines = [line.rstrip('\n') for line in lines]


    if verbose:
        for line in lines[:20]:
            if line[0] == '#' and ":" in line:
                print(line.replace("#", ''))

    if parse_attributes:
        return _parse_all(lines, verbose)

    else:
        return _parse_cols(lines, verbose)