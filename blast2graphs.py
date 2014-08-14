#!/usr/bin/env python

# Standard Python libraries
from sys import stderr
import sys
import argparse
from decimal import Decimal

# Third-party libraries
import networkx as nx


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
    #TODO: validate data columns when headers are present

    met_grf = nx.Graph()  # NetworkX graph with various BLAST-based metrics
    org_ids = set()
    metrics = ['nle', 'bit', 'bsr', 'bal']

    get_self_bit_scores_and_org_ids(met_grf=met_grf, blast_handle=args.blast,
                                    idchar=args.idchar, org_ids=org_ids,
                                    qlcol=args.qlcol-1, slcol=args.slcol-1)

    args.blast.seek(0)

    get_metrics(met_grf=met_grf, blast_handle=args.blast,
                qlcol=args.qlcol-1, slcol=args.slcol-1)

    avgs_wo = compute_organism_averages(
        met_grf=met_grf, metrics=metrics, idchar=args.idchar, org_ids=org_ids)

    compute_global_averages(org_avgs=avgs_wo, metrics=metrics)

    print_unnormalized_abc_files(met_grf=met_grf, metrics=metrics,
                                 glb_avgs=avgs_wo.node['global'],
                                 out_pref=str(args.out_pref)+"_raw")

    print_normalized_abc_files(met_grf=met_grf, metrics=metrics,
                               idchar=args.idchar, org_avgs=avgs_wo,
                               out_pref=str(args.out_pref)+"_nrm")

    if args.fasta:
        print_connected_component_fasta_files(met_grf=met_grf,
                                              fasta_handle=args.fasta,
                                              out_pref=args.out_pref)


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
        description='Generate a set of graphs from a tab-delimited BLASTP ' +
                    'or BLASTN file such that the first two columns contain ' +
                    'the query and subject IDs, respectively, and the last ' +
                    'four columns contain, in order: E-value, bit score, ' +
                    'query length, subject length')

    # Group: IO options
    parser.add_argument('blast', type=argparse.FileType('r'),
                        help='Tab-delimited BLAST file (comment lines are ' +
                             'okay)')
    parser.add_argument('out_pref',
                        help='Prefix for the MCL-compatible "abc" graph files')

    # Group: Formatting options
    parser.add_argument('--qlcol', dest='qlcol',
                        action='store', type=int, default=13,
                        help='One-indexed column containing query lengths ' +
                             '(not required if files include standard ' +
                             'header lines) [def=13]')
    parser.add_argument('--slcol', dest='slcol',
                        action='store', type=int, default=14,
                        help='One-indexed column containing subject lengths ' +
                             '(not required if files include standard ' +
                             'header lines) [def=14]')
    parser.add_argument('--idchar', dest='idchar', action='store', default='|',
                        help='The character used to separate the organism ' +
                             'ID from the rest of the sequence header ' +
                             '[def="|"]')
    parser.add_argument('--fasta', dest='fasta', type=argparse.FileType('r'),
                        help='FASTA file used to generate BLAST results, ' +
                             'will be split into connected components and ' +
                             'reprinted, one file per connected component')

    # Group: TODO
    parser.add_argument('-m', '--merge', dest='merge',
                        action='store_true', default=False,
                        help='Merge sequences from a single organism when ' +
                             'they have non-overlapping alignments to the ' +
                             'same target sequence')

    args = parser.parse_args()

    return args


def get_self_bit_scores_and_org_ids(
        met_grf, blast_handle, idchar=None, org_ids=None,
        evcol=10, bscol=11, qlcol=12, slcol=13):
    """Get bit scores from full-length self-alignments

    Searches an open file for tab-delimited BLAST hit records where the query
    and reference IDs are identical.  It is important that BLAST is run with
    option "--soft_masking true" or the self-alignments are unlikely to be full
    length.

    BLAST does not guarantee particular boundaries for a given alignment, and
    reciprocal alignments frequently have slightly different boundries and
    scores. Neither is more or less valid, so I see no reason to not use the
    larger of the two. I therefore check to make sure subsequent hits do not
    have a greater score than the ones that have already been found. One
    consequence of this is that the intra- and inter-organism average scores
    can not be calculate on the fly.

    Args:
        bsr_graph: A NetworkX graph data structure (does not need to be empty)
        blast_handle: An open file handle containing self-alignments (can
            contain other alignments and/or comment lines beginning with a hash
            '#' character)
        idchar: Character used to delineate between the organism ID and the
            remainder of the sequence ID
        ord_ids: A Python set variable to which organism IDs will be added
        evcol: Column containing BLAST E-values
        bscol: Column containing BLAST bit scores
        qlcol: Column containing query sequence lengths
        slcol: Column containing subject sequence lengths

    Returns:
        Nothing, the NetworkX graph and organsm IDs set data structures are
        edited in place
    """
    for line in blast_handle:
        temp = line.strip().split()
        if not temp:
            continue
        elif temp[0][0] == "#":
            continue

        if temp[0] == temp[1]:
            seq_id = str(temp[0])
            bit_scr = float(temp[bscol])
            org_ids.add(seq_id.split(idchar)[0])

            if not met_grf.has_node(seq_id):
                met_grf.add_node(seq_id, sbs=bit_scr)
            elif bit_scr > met_grf.node[seq_id]['sbs']:
                met_grf.node[seq_id]['sbs'] = bit_scr


def get_metrics(met_grf, blast_handle,
                evcol=10, bscol=11, qlcol=12, slcol=13):
    """Get bit scores from full-length alignments between different sequences

    Searches an open file for tab-delimited BLAST hit records where the query
    and reference IDs are not identical. It is important that BLAST is run with
    option "--soft_masking true" or the self-alignments are unlikely to be full
    length.

    BLAST does not guarantee particular boundaries for a given alignment, and
    reciprocal alignments frequently have slightly different boundries and
    scores. Neither is more or less valid, so I see no reason to not use the
    larger of the two. I therefore check to make sure subsequent hits do not
    have a greater score than the ones that have already been found. One
    consequence of this is that the intra- and inter-organism average scores
    can not be calculated on the fly.

    I convinced myself that removing self-hits from the graph would result in
    more accurate intra-/inter-organism averages and improve normalization.
    They are now removed during printing.

    Args:
        bsr_graph: A NetworkX graph data structure containing self-alignment
            scores
        blast_handle: An open file handle containing non-self-alignments (can
            contain other alignments and/or comment lines beginning with a hash
            '#' character)
        evcol: Column containing BLAST E-values
        bscol: Column containing BLAST bit scores
        qlcol: Column containing query sequence lengths
        slcol: Column containing subject sequence lengths

    Returns:
        Nothing, all data structures are edited in place
    """
    for line in blast_handle:
        temp = line.strip().split()
        if not temp:
            continue
        elif temp[0][0] == "#":
            continue

        metrics = dict()
        qry_id = str(temp[0])
        ref_id = str(temp[1])

        if met_grf.has_node(qry_id) and met_grf.has_node(ref_id):
            qry_len = float(temp[qlcol])
            ref_len = float(temp[slcol])
            aln_len = float(temp[3])
            qry_aln_beg = int(temp[6])
            qry_aln_end = int(temp[7])
            ref_aln_beg = int(temp[8])
            ref_aln_end = int(temp[9])
            metrics['bit'] = float(temp[bscol])

            #BLAST 2.2.28+ rounds E-values smaller than 1e-180 to zero
            if float(temp[evcol]) == 0:
                metrics['nle'] = float(181)
            else:
                # Compute Negative common (base 10) Log of the E-value
                metrics['nle'] = float(-Decimal(temp[evcol]).log10())

            # Compute 'Bit per Anchored length'
            anchored_length = compute_anchored_length(
                qry_aln_beg=qry_aln_beg, qry_aln_end=qry_aln_end,
                ref_aln_beg=ref_aln_beg, ref_aln_end=ref_aln_end,
                aln_len=aln_len, qry_len=qry_len, ref_len=ref_len)
            metrics['bal'] = metrics['bit'] / anchored_length

            # Compute 'Bit Score Ratio'
            qry_sbs = met_grf.node[qry_id]['sbs']
            ref_sbs = met_grf.node[ref_id]['sbs']
            metrics['bsr'] = metrics['bit'] / min(qry_sbs, ref_sbs)

            if not met_grf.has_edge(qry_id, ref_id):
                met_grf.add_edge(qry_id, ref_id)
                for met in metrics.keys():
                    met_grf[qry_id][ref_id][met] = metrics[met]

            # Largest bit score => best hit
            elif metrics['bit'] > met_grf[qry_id][ref_id]['bit']:
                for met in metrics.keys():
                    met_grf[qry_id][ref_id][met] = metrics[met]


def compute_anchored_length(qry_aln_beg, qry_aln_end, ref_aln_beg, ref_aln_end,
                            aln_len, qry_len, ref_len):
    """
    Compute the maximal length of the alignable region anchored by the best hit
    """
    if qry_aln_beg < qry_aln_end:
        qab = qry_aln_beg
        qae = qry_aln_end
    else:
        qab = qry_aln_end
        qae = qry_aln_beg

    if ref_aln_beg < ref_aln_end:
        rab = ref_aln_beg
        rae = ref_aln_end
    else:
        rab = ref_aln_end
        rae = ref_aln_beg

    left_ohang = min(qab, rab)-1
    right_ohang = min(qry_len-qae, ref_len-rae)

    return left_ohang + aln_len + right_ohang


def print_unnormalized_abc_files(met_grf, metrics, glb_avgs, out_pref):
    """Print MCL-formatted .abc graph files"""
    handle = dict()

    for met in metrics:
        handle['dmls_'+met] = open(out_pref+'_dmls_'+met+'.abc', 'w')
        handle['dmnd_'+met] = open(out_pref+'_dmnd_'+met+'.abc', 'w')

    for qry_id, ref_id, edata in met_grf.edges(data=True):
        if qry_id == ref_id:
            continue

        for met in metrics:
            if edata[met]:
                dmls_met = edata[met] / glb_avgs[met+'_avg']
                dmnd_met = edata[met]
            elif met == 'nle':  # Restore p(BLAST-rounded E-value)s
                dmls_met = float(181) / glb_avgs[met+'_avg']
                dmnd_met = float(181)
            else:  # just in case...
                err_out = "Found empty {0} value between {1} and {2}.\n" \
                          .format(met, qry_id, ref_id)
                sys.stderr.write(err_out)

            dmls_out = '{0}\t{1}\t{2}\n'.format(qry_id, ref_id, dmls_met)
            handle['dmls_'+met].write(dmls_out)

            dmnd_out = '{0}\t{1}\t{2}\n'.format(qry_id, ref_id, dmnd_met)
            handle['dmnd_'+met].write(dmnd_out)

    for met in metrics:
        handle['dmls_'+met].close()
        handle['dmnd_'+met].close()


def compute_organism_averages(met_grf, metrics, idchar, org_ids):
    """Compute average scores between and within each pair of organisms

    Args:
        bsr_graph: A NetworkX graph data structure containing a node for each
            sequence and with edges weighted using each BLAST-based metric
        metrics: An ordered list of metrics used in the met_grf data structure
        org_ids: A set containing each organism ID
        idchar: Character used to delineate between the organism ID and the
            remainder of the sequence ID
    Returns:
        avgs_wo: A NetworkX graph data structure containing the total number
            of edges between each pair of organisms, the cumulative sum of each
            metric, and the average score for each metric (one node per
            organism, one edge per pair), excluding (WithOut) self-hits
    """
    avgs_wo = nx.Graph()  # Inter-organism averages w/o self hits

    for qry_id, ref_id, edata in met_grf.edges(data=True):
        qry_org = qry_id.split(idchar)[0]
        ref_org = ref_id.split(idchar)[0]

        # Skip self-hits for the "without" averages
        if not qry_id == ref_id:
            try:
                avgs_wo[qry_org][ref_org]['cnt'] += 1
                for met in metrics:
                    avgs_wo[qry_org][ref_org][met+'_sum'] += edata[met]
            except KeyError:
                avgs_wo.add_edge(qry_org, ref_org, cnt=1)
                for met in metrics:
                    avgs_wo[qry_org][ref_org][met+'_sum'] = edata[met]
                    avgs_wo[qry_org][ref_org][met+'_avg'] = None

    return avgs_wo


def compute_global_averages(org_avgs, metrics):
    """Compute global averages for the entire graph

    Args:
        org_avgs: A NetworkX graph data structure with a node for each organism
            and edges containing the number of edges between each organism in
            the met_grf data structure, and the cumulative weight of all edges
            for each metric between each pair of organisms
        metrics: An ordered list of metrics used in the met_grf data structure
    Returns:
        Nothing, the org_avgs data structure is modified in place
    """
    #FIXME: It would be preferable if the metric_sum names were generated
    #       automatically from the 'metrics' list
    # The 'global' node has degree 0
    org_avgs.add_node('global', cnt=0, bit_sum=float(0), bal_sum=float(0),
                      bsr_sum=float(0), nle_sum=float(0))

    for qry_org, ref_org, edata in org_avgs.edges(data=True):
        org_avgs.node['global']['cnt'] += edata['cnt']

        for met in metrics:
            org_avgs.node['global'][met+'_sum'] += edata[met+'_sum']
            temp_avg = edata[met+'_sum']/edata['cnt']
            org_avgs[qry_org][ref_org][met+'_avg'] = temp_avg

    glb_cnt = org_avgs.node['global']['cnt']

    for met in metrics:
        met_sum = org_avgs.node['global'][met+'_sum']  # floats or Decimals
        org_avgs.node['global'][met+'_avg'] = met_sum/glb_cnt


def print_normalized_abc_files(met_grf, metrics, idchar, org_avgs, out_pref):
    """Normalize metrics to adjust for average inter-organism divergence

    Iterates through the edges in a NetworkX graph, multiplying each metric
    used to weight each edge by a ratio of the average score for that metric
    across the entire graph over the average score for that metric between
    the pair of organisms connected by that particular edge.

    BLAST-rounded E-values and corresponding p(E-values) should be set to zero
    in the metrics graph before calling normalize_metrics(). This function will
    then scale the rest of the graph and identify the largest p(E-value) and
    most minimizing scaling factor. It will then generate a supplemental
    heuristic distribution for the BLAST-rounded values that is guaranteed to
    lie beyond the limits of the real distribution. Lastly, it will convert
    this distribution back into E-values. Supplemental E-values are not
    computed directly because they are more prone to rounding errors and
    because manipulating floats offers significant performance benefits.
    """
    glb_avg = dict()
    handle = dict()

    # Print both DiMensioNeD and DiMensionLesS normalized values
    for met in metrics:
        glb_avg[met] = org_avgs.node['global'][met+'_avg']
        handle['dmls_'+met] = open(out_pref+'_dmls_'+met+'.abc', 'w')
        handle['dmnd_'+met] = open(out_pref+'_dmnd_'+met+'.abc', 'w')

    for qry_id, ref_id, edata in met_grf.edges(data=True):
        if qry_id == ref_id:
            continue

        qry_org = qry_id.split(idchar)[0]
        ref_org = ref_id.split(idchar)[0]

        # I inverted these fractions from what would be more intuitive to avoid
        # some small numbers
        for met in metrics:
            dmls_scl = org_avgs[qry_org][ref_org][met+'_avg']
            dmnd_scl = org_avgs[qry_org][ref_org][met+'_avg'] / glb_avg[met]

            if edata[met]:
                dmls_met = edata[met] / dmls_scl
                dmnd_met = edata[met] / dmnd_scl
            else:  # just in case...
                err_out = "Found empty {0} value between {1} and {2}.\n" \
                          .format(met, qry_id, ref_id)
                sys.stderr.write(err_out)

            dmls_out = '{0}\t{1}\t{2}\n'.format(qry_id, ref_id, dmls_met)
            handle['dmls_'+met].write(dmls_out)

            dmnd_out = '{0}\t{1}\t{2}\n'.format(qry_id, ref_id, dmnd_met)
            handle['dmnd_'+met].write(dmnd_out)

    for met in metrics:
        handle['dmls_'+met].close()
        handle['dmnd_'+met].close()


def print_connected_component_fasta_files(met_grf, fasta_handle, out_pref):
    """
    """
    from Bio import SeqIO
    fasta = SeqIO.to_dict(SeqIO.parse(fasta_handle, 'fasta'))
    w = len(str(len(nx.connected_components(met_grf))))
    cmp_cnt = 0
    for comp in nx.connected_components(met_grf):
        cmp_hdl = open(out_pref+"_comp"+str(cmp_cnt).zfill(w)+".fasta", 'w')
        for sid in comp:
            try:
                cmp_hdl.write(">{0}\n{1}\n".format(sid, fasta[sid].seq))
            except KeyError:
                stderr.write("{0} not found in FASTA file\n".format(sid))
        cmp_hdl.close()
        cmp_cnt += 1


if __name__ == "__main__":
    sys.exit(main())
