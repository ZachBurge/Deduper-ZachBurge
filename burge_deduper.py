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
    fields = line.split()
    # print(fields)
    return (fields[0].split(":")[-1], # qname (specifically umi)
            fields[1], # bitwise flag
            fields[2], # rname
            fields[3], # position
            fields[5]) # cigar

def get_strand(fields):
    flag = int(fields[1])
    if (flag & 16) == 16:
        return "-"
    return "+"

def get_five_prime_pos(fields, strand):
    oblm_pos = int(fields[3]) # 1-based left-most mapping position
    cigar = fields[4]
    parsed_cigar = [(length, op) for length, op in re.findall("(\\d+)([MIDNSHPX=])", cigar)] # split cigar string into tuples of the length and operator for each segment
    if strand == '+':
        leading_clip = 0
        for (length, op) in parsed_cigar:
            if op == 'S':
                leading_clip = int(length)
                break
        return oblm_pos - leading_clip
    else: # strand == '-'
        read_len = 0
        for (length, op) in parsed_cigar:
            if op in ["M", "D", "N", "X", "="]: # operations that consume the reference
                read_len += int(length) # count to the end of the read, since the start will be at the end for the reverse strand
        trailing_clip = 0
        for (length, op) in reversed(parsed_cigar): # reverse the cigar string since the soft clipping would be at the end, a little bit faster
            if op == "S":
                trailing_clip = int(length)
                break
        return oblm_pos + read_len - 1 + trailing_clip

def is_valid_read(fields, known_umis):
    if fields[4] == "*": # no CIGAR, unmapped
        return False
    if fields[0] not in known_umis:
        return False
    return True

def write_pending_reads(pending_reads, outfile):
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
    reads_to_write = []
    seen_keys = set()
    for line in infile:
        line = line.strip()
        if line.startswith("@"):
            outfile.write(f'{line}\n')
        else:
            fields = parse_line(line)

            if fields[2] != curr_chrom:
                write_pending_reads(reads_to_write, outfile)
                seen_keys.clear()
                reads_to_write.clear()
                curr_chrom = fields[2]

            if not is_valid_read(fields, known_umis):
                continue

            umi = fields[0]
            strand = get_strand(fields)
            five_p_pos = get_five_prime_pos(fields, strand)
            key = (fields[2], strand, five_p_pos, umi)

            if key not in seen_keys:
                seen_keys.add(key)
                reads_to_write.append(line)
    if reads_to_write:
        write_pending_reads(reads_to_write, args.outfile)