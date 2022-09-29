#!/usr/bin/env python3
#
'''
This script takes as input a pdb file and marks the requested positions (setting B-factor to 1);
It also uses a fasta with at least one sequence (>$name_long) corresponding to the protein sequence in pdb  
'''

import argparse
import pyfastx
import re
import sys


__author__ = "Ekaterina Osipova, 2022."


def read_fasta_dict(fasta_file):
    ## Reads alignment fasta into a dictionary

    fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fasta_dict[name] = seq
    return fasta_dict


def find_long_short(fasta_dict):
    ## Finds long and short seq

    long_seq = ''
    short_seq = ''
    for k, v in fasta_dict.items():
        if k.endswith('long'):
            long_seq = v
        elif k.endswith('short'):
            short_seq = v
    # check if we got something
    if long_seq == '':
        print('ERROR:Did not find long sequence! at least for one sequence you need to provide _long label!')
        sys.exit(1)
    if short_seq == '':
        print('WARNING:Did not find short sequence: using long sequence as reference')
        short_seq = 'notfound'
    return long_seq, short_seq


def compare_two_seq(seq_long, seq_short):
    ## Makes correspondence between indexes in long and short sequences

    i = 1
    j = 1
    index_to_keep = {}
    for k in range(len(seq_long)):
        if seq_long[k] != '-':
            i += 1
        if seq_short == 'notfound':
            index_to_keep[i] = i
        else:
            if seq_short[k] != '-':
                j += 1
            if (seq_long[k] != '-') and (seq_short[k] != '-'):
                index_to_keep[i] = j
    return index_to_keep


def replace_pdb_bfactor(pdb, sites, index_to_keep, bvalue):
    ## Reads PDB file, replace b-factor for requested sites

    with open(pdb, 'r') as inf:
        for line in inf.readlines():
            if line.startswith('ATOM'):
                res_number = int(line.rstrip()[22:26].rstrip())
                old_pdb_line = line.rstrip()[:60]
                last_field = line.rstrip()[-12:]
                if res_number in index_to_keep:
                    short_index = index_to_keep[res_number]
                    if short_index in sites:
                        new_bfactor = bvalue
                    else:
                        new_bfactor = 0.00
                else:
                    new_bfactor = 0.00
                # output new pdb line
                print('{} {:.2f}{}'.format(old_pdb_line, new_bfactor, last_field))


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='fasta file with protein alignment')
    parser.add_argument('-p', '--pdb', type=str, help='file in pdb format with only ATOM entries')
    parser.add_argument('-s', '--sites', type=str, help='comma-sep residue Nrs to highlight in pdb: 45,278,442')
    parser.add_argument('-b', '--bvalue', type=float, default=1.00,
                        help='new b-factor value you want to highlight the sites with; default=1')
    args = parser.parse_args()

    ## Read fasta
    fasta_dict = read_fasta_dict(args.fasta)

    ## Find long and short seq
    long_seq, short_seq = find_long_short(fasta_dict)

    ## Make site correspondance between long and short sequences
    index_to_keep = compare_two_seq(long_seq, short_seq)

    ## Go through PDB file and assign new b-factor values to corresponding positions in the structure
    site_line = args.sites
    sites = [int(i) for i in  site_line.split(',')]
    replace_pdb_bfactor(args.pdb, sites, index_to_keep, args.bvalue)


if __name__ == "__main__":
        main()