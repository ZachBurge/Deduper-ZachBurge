#!/usr/bin/env python

import argparse
import re

def getargs():
    parser = argparse.ArgumentParser(description="Script to find and remove PCR duplicates from a sorted sam file. \
                                     The script assumes the input sam file is already sorted using samtools sort. \
                                     For any reads that have duplicates, it is the first read that will be kept, and any \
                                     duplicate reads found after are removed. ")
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file")
    parser.add_argument("-o", "--outfile", help="designates absolute file path to output sam file")
    parser.add_argument("-u", "--umi", help="designates file containing the list of UMIs")
    return parser.parse_args()

args = getargs()

def parse_line(line):
    '''
    Function to extract needed information for duplicate checking from sam file line
    '''
    fields = line.split()
    # print(fields)
    return (fields[0].split(":")[-1], # qname (specifically umi)
            fields[1], # bitwise flag
            fields[2], # rname
            fields[3], # position
            fields[5]) # cigar

def get_strand(fields):
    '''
    Function to check the bitwise flag of a sam file line. Checks the 0x16 bit of the flag and returns whether the read is on the '+' or '-' strand.
    '''
    flag = int(fields[1])
    if (flag & 16) == 16:
        return "-"
    return "+"

def get_five_prime_pos(fields, strand):
    '''
    Function to determine the 5' prime position of a read. Uses the 1-based left-most mapping position, the CIGAR string, and the strand of the read to calculate.
    '''
    oblm_pos = int(fields[3]) # 1-based left-most mapping position
    cigar = fields[4]
    parsed_cigar = [(int(length), op) for length, op in re.findall("(\\d+)([MIDNSHPX=])", cigar)] # split cigar string into tuples of the length and operator for each segment
    if strand == '+':
        # 5' end is the left end, so only a soft clip at the beginning matters
        leading_clip = parsed_cigar[0][0] if parsed_cigar[0][1] == "S" else 0
        return oblm_pos - leading_clip
    else: # strand == '-'
        read_len = 0
        for (length, op) in parsed_cigar:
            if op in ["M", "D", "N", "X", "="]: # operations that consume the reference
                read_len += length # count to the end of the read, since the start will be at the end for the reverse strand
        # the 5' end is the right end, so only a soft clip at the end matters
        trailing_clip = parsed_cigar[-1][0] if parsed_cigar[-1][1] == "S" else 0
        return oblm_pos + read_len - 1 + trailing_clip

def is_valid_read(fields, known_umis):
    '''
    Function to check if read is valid to do duplicate checking. Checks for unmapped reads or reads that don't have a known umi.
    '''
    if fields[4] == "*": # no CIGAR, unmapped
        return False
    if fields[0] not in known_umis:
        return False
    return True

def write_pending_reads(pending_reads, outfile):
    '''
    Function to write out pending reads to the outfile. 
    '''
    for read in pending_reads:
        outfile.write(f"{read}\n")

# load known umis from file into set
known_umis = set()
with open(args.umi, "r") as umi_file:
    for line in umi_file:
        line = line.strip()
        known_umis.add(line)

with open(args.file, "r") as infile, open(args.outfile, "w") as outfile:
    curr_chrom = None
    reads_to_write = [] # list of reads that are not duplicates and can be written to the output file
    seen_keys = set() # if a read is valid, build a deduplicate key and add it to this set
    for line in infile:
        line = line.strip()
        if line.startswith("@"): # write out header lines
            outfile.write(f'{line}\n')
        else:
            fields = parse_line(line)

            if fields[2] != curr_chrom: # new chromosome, so write out all previous pending reads and clear everything
                write_pending_reads(reads_to_write, outfile)
                seen_keys.clear()
                reads_to_write.clear()
                curr_chrom = fields[2]

            if not is_valid_read(fields, known_umis): # skip read and go to next
                continue

            umi = fields[0]
            strand = get_strand(fields)
            five_p_pos = get_five_prime_pos(fields, strand)
            key = (fields[2], strand, five_p_pos, umi) # deduplicate key to add to set

            if key not in seen_keys: # unique read, add it to the set and its line to pending reads
                seen_keys.add(key)
                reads_to_write.append(line)
    if reads_to_write: # reads leftover after the loop is done, write them out
        write_pending_reads(reads_to_write, outfile)