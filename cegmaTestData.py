#!/usr/bin/env python

"""
This was forked from the kogTestData.py program.

Rather than make a flexible and reusable Python script, and then write a BASH
wrapper around it to document how I ran it for each data set, I decided to just
hard code this script to generate the data sets I intend to use in my current
study. If someone ever needs to generalize the code, they can grab the useful
functions.
"""

import sys
import argparse

from Bio import SeqIO
from numpy.random import seed, random_integers, shuffle


def main(argv=None):
    """Where the magic happens!"""
    if argv is None:
        argv = sys.argv

    args = get_parsed_args()

    seed(42)

    oe_handle = open("cegma_" + str(args.scheme) + "_ord_evn.fasta", "w")
    or_handle = open("cegma_" + str(args.scheme) + "_ord_rnd.fasta", "w")
    se_handle = open("cegma_" + str(args.scheme) + "_shf_evn.fasta", "w")
    sr_handle = open("cegma_" + str(args.scheme) + "_shf_rnd.fasta", "w")

    cegma = import_cegma(args.fasta)

    for kog_id in sorted(cegma.keys()):

        ord_scheme = list(str(args.scheme))
        shf_scheme = list(str(args.scheme))
        shuffle(shf_scheme)

        for org_id in sorted(cegma[kog_id].keys()):

            pieces = int(ord_scheme.pop())
            even_split(oe_handle, pieces, *cegma[kog_id][org_id])
            rand_split(or_handle, pieces, *cegma[kog_id][org_id])

            pieces = int(shf_scheme.pop())
            even_split(se_handle, pieces, *cegma[kog_id][org_id])
            rand_split(sr_handle, pieces, *cegma[kog_id][org_id])

    oe_handle.close()
    or_handle.close()
    se_handle.close()
    sr_handle.close()


def get_parsed_args():
    """Parse the command line arguments

    Parses command line arguments using the argparse package, which is a
    standard Python module starting with version 2.7.

    Parameters
    ----------
    None, argparse fetches them from user input

    Returns
    -------
    args : argparse.Namespace object
        Contains the parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Split the sequnces in the CEGMA database in order to " +
                    "create several test databases for benchmarking " +
                    "different clustering approaches. Generates four files" +
                    "corresponding to all combinations of two pairs of " +
                    "parameters.  Files with 'even' in the name contain " +
                    "subsequnces of equal size, whereas those with 'rand' " +
                    "contain subsequences of random size.  Files with " +
                    "'ordered' in the name have had the fragmentation " +
                    "scheme applied directly to each organism in " +
                    "alphabetical order (eg. '123456' would leave A. " +
                    "thaliana seqs in tact, split C. elegans seqs in half " +
                    "D. melanogaster seqs into thirds, etc.).  Files with " +
                    "'shuffle' in the name have had the fragmentation " +
                    "scheme shuffled within each KOG in an effort to " +
                    "minimize per-organism effects.")

    parser.add_argument('fasta', type=argparse.FileType('r'),
                        help="Reformatted CEGMA FASTA file")

    parser.add_argument('scheme', type=int,
                        help="A six digit integer encoding the " +
                             "fragmentation scheme usig one digit per " +
                             "organism, where the organism order is: " +
                             "Ath, Cel, ")

    args = parser.parse_args()

    return args


def import_cegma(fasta):
    """Import entire CEGMA database into a Python dictionary

    Python dictionary is keyed by KOG ID at the top level, and by organism ID
    within each KOG.  This organization makes is easier to apply complex
    fragmentation schemes.
    """
    cegma = dict()
    for record in SeqIO.parse(fasta, 'fasta'):
        kog_id, org_id, seq_id, seq_len, seq = seq_info(record)
        try:
            cegma[kog_id][org_id] = (seq_id, seq_len, seq)
        except KeyError:
            cegma[kog_id] = {org_id: (seq_id, seq_len, seq)}

    return cegma


def seq_info(record):
    """Mine info from a SeqIO sequence object

    SeqIO object must have a header with format: orgid|seqid___kogid
    """
    seq_id = str(record.id)
    org_id = seq_id.split('|')[0]
    kog_id = seq_id.split('___')[1]
    seq = str(record.seq)
    seq_len = len(str(record.seq))

    return kog_id, org_id, seq_id, seq_len, seq


def even_split(handle, pieces, seq_id, seq_len, seq):
    """Split sequence into subsequences of equal length"""
    sub_len = int(seq_len/pieces)
    for i in range(pieces):
        header = '>{0}---{1}of{2}'.format(seq_id, i+1, pieces)
        if i+1 < pieces:
            sub_seq = seq[i*sub_len:(i+1)*sub_len]
        else:
            sub_seq = seq[i*sub_len:]
        handle.write("{0}\n{1}\n".format(header, sub_seq))


def rand_split(handle, pieces, seq_id, seq_len, seq):
    """Split sequence into subsequences of random lengths"""
    break_points = get_breaks(seq_len, pieces)
    for i in range(pieces):
        header = '>{0}---{1}of{2}'.format(seq_id, i+1, pieces)
        if i+1 < pieces:
            sub_seq = seq[break_points[i]:break_points[i+1]]
        else:
            sub_seq = seq[break_points[i]:]
        handle.write("{0}\n{1}\n".format(header, sub_seq))


def get_breaks(seq_len, pieces):
    """Select break point positions"""
    internal_breaks = sorted(list(random_integers(2, seq_len, pieces-1)))
    break_points = [0]+internal_breaks+[seq_len]
    bad_break = False
    for i in range(1, len(break_points)):
        if break_points[i]-break_points[i-1] < 3:
            bad_break = True
    if bad_break:
        break_points = get_breaks(seq_len, pieces)

    return break_points


if __name__ == "__main__":
    sys.exit(main())
