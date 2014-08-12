#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 14:29:56 2014

@author: tgibbons
"""

import sys
import argparse
import urllib2
from StringIO import StringIO
import re
from Bio import SeqIO


def main(argv=None):
    """Where the magic happens!

    The main() function coordinates calls to all of the other functions in this
    program in the hope that, by their powers combined, useful work will be
    done.

    I typically try to be pretty good about following PEP-8 guidelines,
    wrapping my code and what not, and not hard-coding information into a
    program without having a runtime option to override it. I work on a pretty
    large screens though, and this program contains many long strings that I
    just found to be much more readable when I let them stretch out.

    The hard-coded information in this program is current as of Aug 7, 2014.

    Args:
        None

    Returns:
        An exit status (hopefully 0)
    """
    if argv is None:
        argv = sys.argv

    three2five = dict(
        ath='Athal',
        cel='Celeg',
        dme='Dmela',
        hsa='Hsapi',
        sce='Scere',
        spo='Spomb',
        ecu='Ecuni')

    # Functional information keyed by single-letter code
    # Adapted from the COG/KOG/fun.txt file
    func = dict(
        J='Information storage and processing: Translation, ribosomal structure and biogenesis',
        A='Information storage and processing: RNA processing and modification',
        K='Information storage and processing: Transcription',
        L='Information storage and processing: Replication, recombination and repair',
        B='Information storage and processing: Chromatin structure and dynamics',
        D='Cellular processes and signaling: Cell cycle control, cell division, chromosome partitioning',
        Y='Cellular processes and signaling: Nuclear structure',
        V='Cellular processes and signaling: Defense mechanisms',
        T='Cellular processes and signaling: Signal transduction mechanisms',
        M='Cellular processes and signaling: Cell wall/membrane/envelope biogenesis',
        N='Cellular processes and signaling: Cell motility',
        Z='Cellular processes and signaling: Cytoskeleton',
        W='Cellular processes and signaling: Extracellular structures',
        U='Cellular processes and signaling: Intracellular trafficking, secretion, and vesicular transport',
        O='Cellular processes and signaling: Posttranslational modification, protein turnover, chaperones',
        C='Metabolism: Energy production and conversion',
        G='Metabolism: Carbohydrate transport and metabolism',
        E='Metabolism: Amino acid transport and metabolism',
        F='Metabolism: Nucleotide transport and metabolism',
        H='Metabolism: Coenzyme transport and metabolism',
        I='Metabolism: Lipid transport and metabolism',
        P='Metabolism: Inorganic ion transport and metabolism',
        Q='Metabolism: Secondary metabolites biosynthesis, transport and catabolism',
        R='Poorly characterized: General function prediction only',
        S='Poorly characterized: Function unknown')

    # Web addresses for the larger files, as of August 8, 2014
    # core.fa is identical to the data/kogs.fa file distributed with CEGMA v2.5
    cegma_url = 'http://korflab.ucdavis.edu/datasets/cegma/core/core.fa'
    kog_url = 'ftp://ftp.ncbi.nih.gov/pub/COG/KOG/kog'
    kyva_url = 'ftp://ftp.ncbi.nih.gov/pub/COG/KOG/kyva'
    Agamb_url = 'http://korflab.ucdavis.edu/datasets/cegma/new_core/A.gambiae.aa'
    Crein_url = 'http://korflab.ucdavis.edu/datasets/cegma/new_core/C.reinhardtii.aa'
    Cinte_url = 'http://korflab.ucdavis.edu/datasets/cegma/new_core/C.intestinalis.aa'
    Tgond_url = 'http://korflab.ucdavis.edu/datasets/cegma/new_core/T.gondii.aa'

    args = get_parsed_args()

    cegma = url2handle(cegma_url)
    ckset = get_cegma_kogs(cegma)
    cegma.close()

    kog = url2handle(kog_url)
    s2k, kog_dat = map_seqs_to_kogs(kog, func, ckset, three2five)
    kog.close()

    kyva = url2handle(kyva_url)
    eck = get_complete_cegma_kogs(kyva, s2k)
    kog.close()

    Agamb = url2handle(Agamb_url)
    add_cegma2_org(eck, "Agamb", Agamb, kog_dat)
    Agamb.close()

    Crein = url2handle(Crein_url)
    add_cegma2_org(eck, "Crein", Crein, kog_dat)
    Crein.close()

    Cinte = url2handle(Cinte_url)
    add_cegma2_org(eck, "Cinte", Cinte, kog_dat)
    Cinte.close()

    Tgond = url2handle(Tgond_url)
    add_cegma2_org(eck, "Tgond", Tgond, kog_dat)
    Tgond.close()

    print_fasta_from_dict(eck, args.out)


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
        description='Generate the Expanded CEGMA KOGs database from the ' +
                    'CEGMA and KOG databases, downloaded from their ' +
                    'respective websites')

    parser.add_argument('out', nargs='?', default='eck.fasta',
                        type=argparse.FileType('w'),
                        help="Optional output filename [def=eck.fasta]")

    args = parser.parse_args()

    return args


def url2handle(url):
    """ Download and open a file

    url2handle downloads a file into a variable and converts it into a file
    handle before returning it.
    """
    try:
        response = urllib2.urlopen(url)
        handle = StringIO(response.read())
    except ValueError:
        sys.stderr("Oops! It looks like {0} has moved or is currently " +
                   "unreachable\n".format(url))
        raise

    return handle


def get_cegma_kogs(cegma):
    """Get list of KOGs used in CEGMA database
    
    get_cegma_kogs iterates through the headers in the FASTA-formatted CEGMA
    database and gathers the set of represented KOGs.
    """
    cks = set()

    for record in SeqIO.parse(cegma, 'fasta'):
        kog = re.search('KOG\d{4}', str(record.id)).group()
        cks.add(kog)

    return cks


def map_seqs_to_kogs(kog, func, ckset, three2five):
    """
    """
    s2k = dict()  # Map sequence IDs to KOG IDs
    kog_dat = dict()

    for line in kog:
        temp = line.strip().split()

        if not temp:
            continue

        elif temp[0][0] == '[':
            kog = str(temp[1])
            if kog in ckset:
                ck = True
                funcs = [' '.join(temp[2:])]
                for f in str(temp[0]):
                    try:
                        funcs.append(func[f])
                    except KeyError:
                        # Skips bracket [] characters
                        pass
    
                description = '; '.join(funcs)
                kog_dat[kog] = description

            else:
                ck = False

        elif ck == True:
            org3 = str(temp[0]).rstrip(':')
            org5 = three2five[org3]
            old_id = str(temp[1])
            new_id = "{0}|{1}___{2} {3}".format(org5, old_id, kog, description)
            s2k[old_id] = new_id

        else:
            pass

    return s2k, kog_dat


def get_complete_cegma_kogs(kyva, s2k):
    """Get complete KOGs from which CEGMA was derived
    
    The CEGMA developers were aiming to identify conserved orthologous
    sequences, so the inparalogs were removed from each candidate KOG cluster.
    Dealing with inparalogs is an important part of homology prediction, 
    however, so we reintroduced them for our study. get_complete_cegma_kogs
    filters the complete KOG database down to just the KOGs used in the
    creation of CEGMA, but retains the inparalogous sequences.
    """
    eck = dict()
    for record in SeqIO.parse(kyva, 'fasta'):
        try:
            newid = s2k[record.id]
            eck[newid] = str(record.seq)
        except KeyError:
            continue

    return eck


def add_cegma2_org(eck, org, c2org, kog_dat):
    """
    """
    for record in SeqIO.parse(c2org, 'fasta'):
        kog = re.search('KOG\d{4}', str(record.id)).group()
        sid = str(record.id).split('|')[0].split('.')[1]
        newid = "{0}|{1}___{2} {3}".format(org, sid, kog, kog_dat[kog])
        eck[newid] = str(record.seq)


def print_fasta_from_dict(dict, out):
    """
    """
    for k, v in dict.iteritems():
        out.write(">{0}\n{1}\n".format(k, v))


if __name__ == "__main__":
    sys.exit(main())
