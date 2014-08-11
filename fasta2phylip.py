#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 20:55:58 2014

@author: tgibbons
"""

import sys
from os import remove
import argparse
from Bio import AlignIO


def main(argv=None):
    """Where the magic happens!

    The main() function coordinates calls to all of the other functions in this
    program in the hope that, by their powers combined, useful work will be
    done.

    Args:
        None

    Returns:
        An exit status (hopefully 0)
    """
    if argv is None:
        argv = sys.argv

    args = get_parsed_args()

    try:
        aln = AlignIO.read(args.fasta, 'fasta')
        AlignIO.write(aln, args.phylip, 'phylip-relaxed')

    except ValueError:
        if args.cleanup:
            fname = args.fasta.name
            pname = args.phylip.name

            args.fasta.close()
            args.phylip.close()

            remove(fname)
            remove(pname)


def get_parsed_args():
    """Parse the command line arguments

    Parses command line arguments using the argparse package, which is a
    standard Python module starting with version 2.7.

    Args:
        None, argparse fetches them from user input

    Returns:
        args: An argparse.Namespace object containing the parsed arguments
    """
    parser = argparse.ArgumentParser(
                 description='Convert FASTA-formatted multiple sequence ' +
                             'alignment into relaxed-Phylip format.')

    parser.add_argument('fasta', type=argparse.FileType('r'),
                        help='FASTA-formatted multiple sequence alignment ' +
                             'input file')
    parser.add_argument('phylip', type=argparse.FileType('w'),
                        help='Phylip-formatted output file')
    parser.add_argument('-c', '--cleanup', action='store_true', default=False,
                        help='Remove empty alignment files')

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    sys.exit(main())
