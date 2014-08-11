#!/usr/bin/env python


import sys
import argparse
from copy import deepcopy

from Bio import SeqIO
from numpy.random import seed, random_integers, shuffle


def main(argv=None):
    """Where the magic happens!

    The main() function coordinates calls to all of the other functions in this
    program in the hope that, by their powers combined, useful work will be
    done.
    """
    if argv is None:
        argv = sys.argv

    args = get_parsed_args()

    seed(42)
    eck, orgs, min_len = import_fasta(args.fasta)
    scheme = fit_string_to_length(str(args.scheme), len(orgs))
    scheme_list = [int(x) for x in scheme]
    min_frag = calculate_minimum_fragment_length(min_len, max(scheme_list))

    handle_pref = str(args.prefix) + '_' + scheme
    oe_handle = open(handle_pref + "_ord_evn.fasta", "w")
    or_handle = open(handle_pref + "_ord_rnd.fasta", "w")
    se_handle = open(handle_pref + "_shf_evn.fasta", "w")
    sr_handle = open(handle_pref + "_shf_rnd.fasta", "w")

    for kog_id in sorted(eck.keys()):

        ord_scheme = deepcopy(scheme_list)
        shf_scheme = shuffled_scheme(scheme_list, orgs, eck[kog_id])

        for org_id in sorted(eck[kog_id].keys()):
            ord_pieces = ord_scheme.pop()

            for seq_id in sorted(eck[kog_id][org_id].keys()):
                # seq_info contains a tuple that is unpacked by the * operator
                seq_info = eck[kog_id][org_id][seq_id]

                even_split(oe_handle, ord_pieces, seq_id, *seq_info)
                rand_split(or_handle, ord_pieces, min_frag, seq_id, *seq_info)

                shf_pieces = shf_scheme.pop()
                even_split(se_handle, shf_pieces, seq_id, *seq_info)
                rand_split(sr_handle, shf_pieces, min_frag, seq_id, *seq_info)

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
        description="Split the sequences from any number of variants of the " +
                    "COG database in order to create several test databases " +
                    "for benchmarking different clustering approaches. " +
                    "Generates four files corresponding to all combinations " +
                    "of two pairs of parameters.  Files with 'even' in the " +
                    "name contain subsequnces of equal size, whereas those " +
                    "with 'rand' contain subsequences of random size.  " +
                    "Files with 'ordered' in the name have had the " +
                    "fragmentation scheme applied directly to each organism " +
                    "in alphabetical order (eg. '123456' would leave A. " +
                    "thaliana seqs in tact, split C. elegans seqs in half " +
                    "D. melanogaster seqs into thirds, etc.).  Files with " +
                    "'shuffle' in the name have had the fragmentation " +
                    "scheme shuffled within each KOG in an effort to " +
                    "minimize per-organism effects.")

    parser.add_argument('fasta', type=argparse.FileType('r'),
                        help="Reformatted COG, KOG, or CEGMA FASTA file")

    parser.add_argument('scheme', type=int,
                        help="An integer encoding the fragmentation scheme " +
                             "usig one digit per organism, where the " +
                             "organism order is: Agamb, Athal, Celeg, " +
                             "Cinte, Crein, Dmela, Ecuni, Hsapi, Scere, " +
                             "Spomb, Tgond. Integers that are not eleven " +
                             "digits will be either repeated or truncated, " +
                             "as necessary.")

    parser.add_argument('prefix', default="eck",
                        help="Prefix for output files [def='eck']")

    args = parser.parse_args()

    return args


def import_fasta(fasta):
    """Import entire ECK database into a Python dictionary

    Python dictionary is keyed by KOG ID at the top level, and by organism ID
    within each KOG.  This organization makes is easier to apply complex
    fragmentation schemes.
    """
    eck = dict()
    orgs = set()
    min_len = float("inf")

    for record in SeqIO.parse(fasta, 'fasta'):
        kog_id, org_id, seq_id, seq_len, seq = sequence_info(record)
        orgs.add(org_id)

        try:
            eck[kog_id][org_id][seq_id] = (seq_len, seq)
        except KeyError:
            try:
                eck[kog_id][org_id] = {seq_id: (seq_len, seq)}
            except KeyError:
                eck[kog_id] = {org_id: {seq_id: (seq_len, seq)}}

        if seq_len < min_len:
            min_len = seq_len

    orgs = list(sorted(orgs))

    return eck, orgs, min_len


def sequence_info(record):
    """Mine info from a SeqIO sequence object

    SeqIO object must have a header with format: orgid|seqid___kogid
    """
    seq_id = str(record.id)
    org_id = seq_id.split('|')[0]
    kog_id = seq_id.split('___')[1]
    seq = str(record.seq)
    seq_len = len(str(record.seq))

    return kog_id, org_id, seq_id, seq_len, seq


def fit_string_to_length(string, length):
    """Repeat string up to specified length

    Modified from a post by Jason Scheirer:
    http://stackoverflow.com/questions/3391076/repeat-string-to-certain-length
    """
    return (string*((length/len(string))+1))[:length]


def calculate_minimum_fragment_length(min_len, max_frags):
    """
    """
    min_even = min_len / max_frags
    return min(10, min_even)


def shuffled_scheme(scheme_list, orgs, kog_dict):
    """Determine a per-COG scheme for fragmenting sequences

    In all but the CEGMA set, the number of represented organisms and the
    number of sequences contributed by each organism varies between COG
    clusters.  This function matches the user-defined fragmentation scheme to
    the organisms represented within a particular COG, counts the number of
    sequences contributed to that COG cluster by that organism, and adds a
    particular fragmentation count to the shuffled scheme as many times as
    there were sequences contributed from the corresponding organism.  It then
    shuffles the order of fragmentation counts before returning the list.
    """
    shf_scheme = list()
    for i in range(len(orgs)):
        frags = [scheme_list[i]]
        try:
            shf_scheme.extend(frags * len(kog_dict[orgs[i]].keys()))
        except KeyError:
            continue

    shuffle(shf_scheme)

    return shf_scheme


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


def rand_split(handle, pieces, min_frag, seq_id, seq_len, seq):
    """Split sequence into subsequences of random lengths"""
    break_points = get_breaks(min_frag, seq_len, pieces)
    for i in range(pieces):
        header = '>{0}---{1}of{2}'.format(seq_id, i+1, pieces)
        if i+1 < pieces:
            sub_seq = seq[break_points[i]:break_points[i+1]]
        else:
            sub_seq = seq[break_points[i]:]
        handle.write("{0}\n{1}\n".format(header, sub_seq))


def get_breaks(min_frag, seq_len, pieces):
    """Select break point positions"""
    internal_breaks = sorted(list(random_integers(2, seq_len, pieces-1)))
    break_points = [0]+internal_breaks+[seq_len]
    bad_break = False
    for i in range(1, len(break_points)):
        if break_points[i]-break_points[i-1] < min_frag:
            bad_break = True
    if bad_break:
        break_points = get_breaks(min_frag, seq_len, pieces)

    return break_points


if __name__ == "__main__":
    sys.exit(main())

