#!/usr/bin/env python3

import argparse
from argparse import RawTextHelpFormatter
import re
import sys
import logging

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description="Subsets an input FASTA file (i.e., a whole genome) based on the output of a BLAST search. Returns a FASTA file containing all of the BLAST hits. Optionally filters the BLAST hits based on arbitrary criteria.", formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--fasta", help="String. Path to input FASTA file containing a genome's contigs.\nDefault: None", required=True)
parser.add_argument("-b", "--blastoutput", help="String. Path to BLAST output table. Expects a BLAST output table specified with the parameters:\n-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovs evalue bitscore\"\nDefault: None", required=True)
parser.add_argument("-f", "--filter", help="String. Filters applied to BLAST table, separated by semicolon and wrapped in quotes. Filters are specified using the format \"'attribute''equals/greater than/less than''value'\"\nFor example: \"qseqid=SpecificAccession; pident>0.90; bitscore>300\"\nIf \"=\" is the comparison, it does simple exact string matching.\nDefault: None", required=False)
parser.add_argument("-o", "--output", help="String. Path to store the regions extracted from the input genome matching the BLAST table after filtering.\nDefault: stdout", required=False)
parser.add_argument("-filtered_table", "--filtered_table", help="String. Path to store the filtered BLAST table in. Mostly useful for debugging purposes.\nDefault: None", required=False)
args = parser.parse_args()

def read_fasta_file(fasta_file):
    """
    Reads a fasta file into a dictionary. Keys are FASTA headers, values are sequences (newlines removed).
    """
    fasta_dictionary = {}
    with open(fasta_file, "r") as fasta_handler:
        for line in fasta_handler:
            line = line.strip()
            if line.startswith(">"):
                header = line.split(">")[1]
                fasta_dictionary[header] = []
            else:
                fasta_dictionary[header].append(line)
    fasta_dictionary = {header : "".join(fasta_dictionary[header]) for header in fasta_dictionary}
    return(fasta_dictionary)

def read_blast_table(blast_file):
    """
    Read in a BLAST hit file. Keys are line numbers, values are BLAST hit records.
    """
    entry = 0
    blast_dict = {}
    with open(blast_file) as blast_table_handler:
        for line in blast_table_handler:
            line = line.strip()
            entry += 1
            blast_dict[entry] = { blast_attributes[i] : line.split("\t")[i] for i in range(len(blast_attributes))}
    return(blast_dict)

def evaluate_filter(blast_entry, filter_statement):
    """
    Looks at a line in a BLAST output table and evaluates an arbitrary filter statement.
    """
    filter_attribute = filter_statement[0]
    filter_comparison = filter_statement[1]
    filter_critical_value = filter_statement[2]

    blast_entry_value = blast_entry[filter_attribute]

    match_status = False
    if filter_comparison == "=":
        if str(filter_critical_value) == str(blast_entry_value):
            match_status = True
    elif filter_comparison == ">":
        if float(blast_entry_value) > float(filter_critical_value):
            match_status = True
    elif filter_comparison == "<":
        if float(blast_entry_value) < float(filter_critical_value):
            match_status = True
    return match_status

def filter_blast_table(blast_table_dict, filter_list):
    """
    Takes a BLAST hit dictionary and returns a dictionary where every record passes some
    arbitrary filter criteria.
    """
    filtered_blast_table_dict = {}

    for key in blast_table_dict:
        blast_entry = blast_table_dict[key]
        all_match = True
        for filter_statement in filter_list:
            check_filter = evaluate_filter(blast_entry, filter_statement)
            if not check_filter:
                all_match = False
                break
        if all_match:
            #print(blast_entry, filter_list)
            filtered_blast_table_dict[key] = blast_entry
    return filtered_blast_table_dict

def reverse_complement(sequence):
    """
    Reverse-complements a DNA sequence.
    """
    rc_dict = {"a":"t", "t":"a", "c":"g", "g":"c", "A":"T", "T":"A", "C":"G", "G":"C", "n":"n", "N":"N"}
    reverse_complement = "".join(rc_dict.get(base, base) for base in reversed(sequence))
    return reverse_complement

def retrieve_subsequences(fasta_dict, filtered_blast_table_dict):
    """
    Given a dictionary of a genome assembly and a list of BLAST hits, extracts the relevant hit regions,
    reverse complementing them as necessary.
    """
    filtered_fasta_dict = {}

    for key in filtered_blast_table_dict:
        blast_hit = filtered_blast_table_dict[key]
        sseqid = blast_hit["sseqid"]
        sstart = blast_hit["sstart"]
        send = blast_hit["send"]
        for fasta_header in fasta_dict:
            fasta_sequence = fasta_dict[fasta_header]
            simple_fasta_header = fasta_header.split(" ")[0]
            if sseqid == simple_fasta_header:
                if sstart < send:
                    # BLAST counts from 1, not 0. So, substract 1 from the indices.
                    subsequence = fasta_sequence[int(sstart)-1 : int(send)-1]
                else:
                    subsequence = reverse_complement(fasta_sequence[int(send)-1 : int(sstart) - 1])
                mod_fasta_header = fasta_header + " " + sstart + "-" + send
                filtered_fasta_dict[mod_fasta_header] = subsequence
    return filtered_fasta_dict

def write_filtered_table(filtered_blast_table_dict, filtered_blast_table_file=False):
    """
    Function for writing an output table of sequences that pass the filter provided to the script.
    Mostly useful for debugging purposes.
    """
    if filtered_blast_table_file:
        with open(filtered_blast_table_file, "w") as blast_handler:
            for key in filtered_blast_table_dict:
                output_line = "\t".join(filtered_blast_table_dict[key].values())
                blast_handler.write("%s\n" % output_line)
    else:
        pass

def write_filtered_fasta(filtered_fasta_dict, filtered_fasta_file=False):
    """
    Function for writing the subsequences from the BLAST hit table.
    """
    if filtered_fasta_file:
        with open(filtered_fasta_file, "w") as fasta_handler:
            for fasta_header in filtered_fasta_dict:
                fasta_handler.write(">%s\n%s\n" % (fasta_header, filtered_fasta_dict[fasta_header]))
    else:
        for fasta_header in filtered_fasta_dict:
            print(">" + fasta_header)
            print(filtered_fasta_dict[fasta_header])

def subset_fasta_by_BLAST_results(filtered_fasta_dict, filtered_blast_table_dict, filtered_fasta_file=False, filtered_blast_table_file=False):
    """
    A wrapper function to call the two writing functions.
    """
    write_filtered_fasta(filtered_fasta_dict, filtered_fasta_file=filtered_fasta)
    write_filtered_table(filtered_blast_table_dict, filtered_blast_table_file=filtered_table)

if args.fasta:
    fasta_file = args.fasta
if args.blastoutput:
    blast_file = args.blastoutput
filtered_table = False
if args.filtered_table:
    filtered_table = args.filtered_table
filtered_fasta = False
if args.output:
    filtered_fasta = args.output

blast_attributes = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send","qcovs", "evalue", "bitscore"]
filter_list = []

if args.filter:
    # remove spaces; split into groups separated by semicolon
    all_filters = args.filter.replace(" ", "").split(";")
    # remove empty instances so that it can be split with a regular expression easily
    all_filters = list(filter(None, all_filters))
    comparison_symbols = re.compile('[>=<]')
    for f in all_filters:
        try:
            attribute = re.split(comparison_symbols, f)[0]
            comparison = re.findall(comparison_symbols, f)[0]
            critical_value = re.split(comparison_symbols, f)[1]
            # Make sure the attribute specified actually exists
            if attribute not in blast_attributes:
                raise
            # Make sure that the argument after ">" or "<" is actually a number.
            if comparison in [">", "<"]:
                try:
                    float(critical_value)
                except:
                    raise
            # Make sure that attributes are not set more than once
            if attribute not in [a[0] for a in filter_list]:
                filter_list.append((attribute, comparison, critical_value))
            else:
                raise
        except:
            # Print a non-specific help message telling the user their filter options were not parsed correctly.
            # Then catch on fire.
            print("\nSomehow your filter cannot be evaluated.\nMake sure you are specifying it correctly (try using --help),\nand make sure to not specify the same option more than once.\n\nYour filter options:\n\t%s\nThe list of acceptable attributes includes:\n\t%s\n" % (all_filters, blast_attributes))
            sys.exit(1)
else:
    # My quick recommendations: ("qcovs", ">", 70), ("evalue", "<", 1e-5)
    # That filter quickly identifies good matches (by E-value) in length-appropriate queries
    # (by using the qcovs value).
    filter_list = []

fasta_dict = read_fasta_file(fasta_file)
blast_table_dict = read_blast_table(blast_file)
filtered_blast_table_dict = filter_blast_table(blast_table_dict, filter_list)

if len(filtered_blast_table_dict) > 0:
    filtered_fasta_dict = retrieve_subsequences(fasta_dict, filtered_blast_table_dict)
else:
    # If there are no entries after filtering, throw an exception.
    print("Error: After filtering, there are no sequences left! Double check your filter statement.")
    sys.exit(1)

if __name__ == "__main__":
    subset_fasta_by_BLAST_results(filtered_fasta_dict, filtered_blast_table_dict, filtered_fasta_file=filtered_fasta, filtered_blast_table_file=filtered_table)
