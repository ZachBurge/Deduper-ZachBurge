Given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory. 

**Define the problem:**  PCR duplicates are identical strands of DNA/RNA that arise from doing PCR amplification during sequencing. They are important to remove because they could result in an incorrect measure of sequencing depth for certain strands of DNA/RNA over others, and many analytics that result from this sequencing data rely on the counts of each strand after sequencing. The goal of this project is to identify and remove PCR duplicates from a SAM file of uniquely mapped reads. 

PCR duplicates will:

* Be on the same chromosome (Field 3 in SAM header)
* Have the same left-most mapping position (Field 4 in SAM header) and same 5' position (to account for soft clipping)
* Be on the same strand (+ or -)
* Have the same UMI (QNAME, SAM col 1)

**Algorithm:**

1) Run samtools sort on the file
2) Load in the UMIs from the file into some object called umis (maybe a list, maybe a set, not sure which is better)
3) Open input sam file and output sam file
4) Initialize pending_reads, a list to contain all of the valid reads ready to be written to the output file.
Also initialize seen_keys (set) to track record information and ensure only 1 of each read makes it to the output file.
5) current_chrom = NULL
6) For line in file:
7) If line starts with '@', it's a header and can just be written to the output file
8) Read in sam line and parse it, called record (parse_line)
9) if rname == current_chrom, write out all records in pending_reads and then empty pending_reads. current_chrom = rname
10) if not is_valid_record(record, umis), skip that record
11) Create the key for the record: umi = get_umi(qname), strand = get_strand(flag), five_prime_pos = get_five_prime_pos(position, cigar, strand),
key = (rname, strand, five_prime_pos, umi)
12) If that key has not been seen (not in seen_keys), add the key to seen keys and append the line to pending_reads. 
13) After the loop, check if there are reads in pending_reads and write them out if there are

Functions:

```
parse_line(line):
'''Function to read in a sam line and extract the needed information for duplicate checking later on'''
fields = line.strip().split()
return (qname(fields[0]), flag(fields[1]), rname(fields[2]), position(fields[3]), cigar(fields[5]))

is_valid_record(record, umi):
'''Function to check if a record is a PCR duplicate'''
Check flag for reverse or not, check umi to see if it is in the known umis
return true if all are true, false if any are false

get_umi(qname):
return last element of split qname on ':'

get_strand(flag):
if flag & 16 == 16, return '-', else return '+'

get_five_prime_pos(position, cigar, strand):
split cigar string into each section (i.e. 2S12M) -> (2, 'S'), (12, 'M')
for each tuple in split cigar string:
    if tuple[1] == 'S':
        soft_clipped = tuple[0]
return position - soft_clipped

```