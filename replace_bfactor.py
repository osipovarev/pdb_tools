#!/usr/bin/env python3
#
import argparse

__author__ = "Ekaterina Osipova, 2020."


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pdb', type=str, help='file in pdb format with only ATOM entries')
    parser.add_argument('-r', '--replace', type=str, help='file with a list of values to replace 11th column in pdb')
    #parser.add_argument('-ss', '--startstop', type=str, defaulthelp='')
    args = parser.parse_args()

    ## Read replacement file into a list
    replace_values = []
    with open(args.replace, 'r') as inf:
        for line in inf.readlines():
            replace_values.append(line.rstrip()[:4])

    ## Read PDB file, replace 11th column with values from the list
    with open(args.pdb, 'r') as inf:
        first_line = inf.readline()
        # -c24-26
        res_number = first_line.rstrip()[23:26]
        i = 0
        old_pdb_line = first_line.rstrip()[:61]
        last_field = first_line.rstrip()[-12:]
        new_bfactor = replace_values[i]
        # output new pdb line
        print('{}{}{}'.format(old_pdb_line, new_bfactor, last_field))

        for line in inf.readlines():
            if line.startswith('ATOM'):
                new_res_number = line.rstrip()[23:26]

                if new_res_number != res_number:
                    i += 1
                    res_number = new_res_number
                old_pdb_line = line.rstrip()[:61]
                last_field = line.rstrip()[-12:]
                new_bfactor = replace_values[i]
                # output new pdb line
                print('{}{}{}'.format(old_pdb_line, new_bfactor, last_field))


if __name__ == "__main__":
    main()