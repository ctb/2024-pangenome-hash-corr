#! /usr/bin/env python
"""
A hackity-hack script used to run through a gather file and extract only
those names that have a certain abundance match in them.

CTB 5/31/2024
"""
import sys
import csv
import argparse
import os.path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gather_csvs', nargs='+')
    p.add_argument('--min-abund', type=float, default=5)
    args = p.parse_args()

    keeplist = []
    for gather_csv in args.gather_csvs:
        with open(gather_csv, 'r', newline='') as fp:
            for row in csv.DictReader(fp):
                med_abund = row['median_abund']
                if med_abund:
                    med_abund = float(med_abund)
                    if med_abund >= args.min_abund:
                        basename = os.path.basename(gather_csv)
                        keeplist.append(basename.split('.')[0])

    for k in keeplist:
        print(k)
    

if __name__ == '__main__':
    sys.exit(main())
