#!/usr/bin/env python3

import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import warnings

def parse_args(args=None):
    Description = "Re-labels the headers of an input FASTA file. Because long headers sometimes become non-unique when transitioning between file formats (e.g., FASTA to Phylip), this simple tool renames headers to simple integers. A TSV table is generated, which can be used to regenerate the original the FASTA file."
    Epilog = "Example usage: python3 fasta_headers_to_numbers.py --fasta [fasta file] --renamed-fasta [renamed fasta file] --renamed-tsv [renamed fasta table]"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-f", "--fasta", default=False, required=False, help="String. Path to input FASTA file. Default: None")
    parser.add_argument("-r", "--renamed-fasta", default=None, required=True, help="String. Path to the renamed output FASTA file. Default: None")
    parser.add_argument("-t", "--renamed-tsv", default=None, required=True, help="String. Path to the renamed output TSV file. Default: None")
    parser.add_argument("--regenerate", default=False, action="store_true", help="Boolean. Operates the script in reverse: instead of enumerating and renaming the headers in an initial FASTA file, use the outputs of a previous iteration of the script (renamed FASTA and TSV files) to regenerate the original FASTA file.")
    return parser.parse_args(args)

def read_fasta_file(fasta_file):
    """
    Reads a fasta file into a dictionary.

    Note that the structure of the data here differs from the output dictionary structure. Keys here are sequence counts, and the values are lists with the elements [fasta header, sequence].
    """

    seq_count = 0
    fasta_dictionary = {}
    with open(fasta_file, "r") as fasta_handler:
        for line in fasta_handler:
            line = line.strip()
            if line.startswith(">"):
                seq_count += 1
                header = line.split(">")[1]
                fasta_dictionary[seq_count] = [header, ""]
            else:
                fasta_dictionary[seq_count][1] = fasta_dictionary[seq_count][1] + line
    return(fasta_dictionary)

def input_dict_to_simple_dict(input_dict):
    """
    A function to turn an input dictionary (i.e., {seq_num : [header, seq] }) into a simplified dictionary (i.e., {header : seq}).
    """
    simple_dict = {}
    for seq_num in input_dict:
        original_header = input_dict[seq_num][0]
        seq = input_dict[seq_num][1]
        if original_header not in simple_dict:
            simple_dict[original_header] = seq
        else:
            warnings.warn("Warning! The FASTA headers in the input FASTA file do not seem to be unique. It is likely that this script will produce incorrect or incomplete information!")
    return simple_dict

def regenerate_original_fasta(renamed_tsv_file, renamed_fasta_file):
    """
    This function uses renamed TSV and renamed FASTA files (i.e., the output of this script) to re-generate the original FASTA input.

    Only needed when --reverse is True.
    """
    tsv_dict = {}
    with open(renamed_tsv_file) as tsv_handler:
        for line in tsv_handler:
            line = line.strip()
            tsv_dict[line.split("\t")[0]] = line.split("\t")[1]
    fasta_dict = read_fasta_file(renamed_fasta_file)

    regenerated_dict = {fasta_dict[seq_count][0] : fasta_dict[seq_count][1] for seq_count in tsv_dict}

    return regenerated_dict

def write_fasta(simple_dict, renamed_fasta_file):
    """
    A simple function to write a FASTA file. A "simple_dict" is a dictionary with {fasta_header : sequence} key/value pairs.
    """
    with open(renamed_fasta_file, "w") as fasta_handler:
        for header in simple_dict:
            fasta_handler.write(">%s\n%s\n" % (header, simple_dict[header]))

def write_renamed_tsv(input_dict, renamed_tsv_file):
    """
    Turns an input FASTA into a TSV with the columns "seq_count[tab]original_header".
    """
    with open(renamed_tsv_file, "w") as tsv_handler:
        for seq in input_dict:
            tsv_handler.write("%s\t%s\n" % (seq, input_dict[seq][0]))

def simplify_fasta_headers(fasta, renamed_fasta, renamed_tsv, regenerate=False):
    if fasta:
        input_dict = read_fasta_file(fasta)
        simple_dict = input_dict_to_simple_dict(input_dict)
        write_renamed_tsv(input_dict, renamed_tsv_file=renamed_fasta)
        write_fasta(simple_dict, renamed_fasta_file=renamed_tsv)
    else:
        if regenerate:
            simple_dict = regenerate_original_fasta(renamed_fasta_file=renamed_fasta, renamed_tsv_file=renamed_tsv)
            write_fasta(simple_dict, renamed_fasta_file="orig."+os.path.basename(renamed_fasta))
        else:
            warnings.warn("Warning! Renamed FASTA and TSV files were provided, but --regenerate was not specified.")

def main(args=None):
    args = parse_args(args)
    simplify_fasta_headers(fasta=args.fasta, renamed_fasta=args.renamed_fasta, renamed_tsv=args.renamed_tsv, regenerate=args.regenerate)

if __name__ == "__main__":
    sys.exit(main())
