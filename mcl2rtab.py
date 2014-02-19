#!/usr/bin/env python

"""
This script should work with clusters of any sequences whose IDs contain the
typical KOG identifiers in their standard format, the keyword 'KOG' followed
by a four digit ID number (eg. KOG0001, KOG2437, etc.).
"""

import sys
import re
import argparse
#from pandas import DataFrame


def main(argv=None):
    """Where the magic happens!

    The main() function coordinates calls to all of the other functions in this
    program in the hope that, by their powers combined, useful work will be
    done.
    """
    if argv is None:
        argv = sys.argv

    args = get_parsed_args()

    kpc_handle = open(args.prefix+"_kogs_per_cluster_summary.Rtab", 'w')
    kpc_handle.write("Order\tFragmentation\tEvalueCutoff\tNormalization\t" +
                     "Metric\tInflation\tKOGsPerCluster\tClusterCount\n")

    cpk_handle = open(args.prefix+"_clusters_per_kog_summary.Rtab", 'w')
    cpk_handle.write("Order\tFragmentation\tEvalueCutoff\tNormalization\t" +
                     "Metric\tInflation\tClustersPerKOG\tKOGCount\n")

#    cdf_cols = ["Metric", "Normalized", "Inflation", "KOGsPerCluster",
#                "ClusterCount"]
#    cluster_df = DataFrame(columns=cdf_cols)
#
#    kdf_cols = ["Metric", "Normalized", "Inflation", "ClustersPerKOG",
#                "KOGCount"]
#    kog_df = DataFrame(columns=kdf_cols)

    for mcl_file in args.mcl_files:
        ordr, frag, ctof, norm, mtrc, infl = parse_file_name(mcl_file.name)
        kogs_per_cluster, clusters_per_kog = score_clustering(mcl_file)

        print_kpc(
            kpc_handle, ordr, frag, ctof, norm, mtrc, infl, kogs_per_cluster)

        print_cpk(
            cpk_handle, ordr, frag, ctof, norm, mtrc, infl, clusters_per_kog)

#        cluster_df(ordr, frag, ctof, norm, mtrc, infl, kogs_per_cluster)
#        kog_df(ordr, frag, ctof, norm, mtrc, infl, clusters_per_kog)

    kpc_handle.close()
    cpk_handle.close()

#    print_cluster_summary(args.prefix, cluster_df)
#    print_kog_summary(args.prefix, kog_df)


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
        description="Generate a set of summary statistics for a given set " +
                    "of MCL clustering files")

    parser.add_argument('prefix',
                        help='Prefix for global summary files')
    parser.add_argument('mcl_files', nargs='+', type=argparse.FileType('r'),
                        help='MCL output files')

    args = parser.parse_args()

    return args


def parse_file_name(mcl_file_name):
    """Gather information about cluster generation parameters from file name
    """
    # Determine of fragmentation scheme was randomized for each KOG
    if re.search('_ord', mcl_file_name):
        ordr = "Ordered"
    elif re.search('_shf', mcl_file_name):
        ordr = "Shuffled"
    else:
        raise Exception(
            "Could not determine if fragmentation scheme was shuffled or " +
            "applied directly for file "+mcl_file_name+". Make sure file " +
            "names contain either '_ord' or '_shf'.")

    # Determine if break points were evenly spaced
    if re.search('_evn', mcl_file_name):
        frag = "Even"
    elif re.search('_rnd', mcl_file_name):
        frag = "Random"
    else:
        raise Exception(
            "Could not determine if break points were evenly spaced in " +
            "file "+mcl_file_name+". Make sure file names contain either " +
            "'_evn' or '_rnd'.")

    # Determine E-value cutoff used for BLASTp
    evl_ctof = re.search('_1e-\d+', mcl_file_name)
    if evl_ctof:
        ctof = evl_ctof.group(0).lstrip('_')
    else:
        raise Exception(
            "Could not determine BLAST E-value cutoff for file " +
            mcl_file_name+". Make sure file names indicate the E-value " +
            "cutoff using fomat '_1e-X'.")

    # Determine if graph was normalized (by inter-/intr-organsm averages)
    if re.search('_nrm', mcl_file_name):
        norm = "Normalized"
    elif re.search('_raw', mcl_file_name):
        norm = "Raw"
    else:
        raise Exception(
            "Could not determine if file "+mcl_file_name+" was normalized. " +
            "Make sure file names contain either '_nrm' or '_raw'.")

    # Identify metric used to weight graph
    if re.search('_bit', mcl_file_name):
        mtrc = "BitScore"
    elif re.search('_bpr', mcl_file_name):
        mtrc = "BitPerResidue"
    elif re.search('_bsr', mcl_file_name):
        mtrc = "BitScoreRatio"
    elif re.search('_evl', mcl_file_name):
        mtrc = "Evalue"
    else:
        raise Exception(
            "Could not determine metric used for file "+mcl_file_name+". " +
            "Make sure file names contain one of '_bit', '_bpr', '_bsr', " +
            "or '_evl'.")

    # Identify inflation parameter used by MCL
    if re.search('I\d{2}', mcl_file_name):
        infl = re.search('_I\d{2}', mcl_file_name).group()
        infl = float(infl.lstrip('_I'))/10
    else:
        raise Exception(
            "Could not determine inflation value for "+mcl_file_name+". " +
            "Make sure file names contain '_I##' where '##' is the " +
            "inflation parameter (sans decimal)")

    return ordr, frag, ctof, norm, mtrc, infl


def score_clustering(mcl_file):
    """Gather statistics from an MCL cluster file

    The function focuses on two measures of success:

    1) KOGs per cluster: Ideally, each cluster will contain sequences from only
        a single KOG. When this is not true, the clustering has combined
        sequences that human experts believe should be clustered separately.
        This measure is analogous to specificity.

    2) Clusters per KOG: Ideally, all members of a KOG will be contained within
        a cluster. When this is not true, the clustering has separted
        sequences that human experts believe should be combined. This measure
        is analogous to sensitivity.

    Parameters
    ----------
    mcl_file : readable_file_handle
        A set of MCL clusters

    Returns
    -------
    kogs_per_cluster : dict
        Keyed by the number of KOGs represented within a particular cluster,
        and storing the number of clusters containing members of that number of
        KOGs as values.
    clusters_per_kog : dict
        Keyed by the number of clusters containing members of a particular KOG,
        and storing the number of KOGs spread across that number of clusters as
        values.
    """
    cluster_counts_by_kog = dict()  # clusters/KOG requires 2-step processing
    clusters_per_kog = dict()
    kogs_per_cluster = dict()

    # Print a per-line KOG summary to a special file for each MCL clustering
    per_cluster_stats = open(mcl_file.name+"-kog_summary", 'w')

    # Parse each cluster
    for cluster in mcl_file:
        seqs = cluster.strip().split()
        kog_counts = dict()

        # Count occurance of each KOG within cluster
        for seq in seqs:
            kog = re.search('KOG\d{4}', seq).group()
            try:
                kog_counts[kog] += 1
            except KeyError:
                kog_counts[kog] = 1

        # Print occurance of each KOG within cluster
        out_line = ''
        for k, v in sorted(kog_counts.iteritems()):
            out_line += str(k)+':'+str(v)+'\t'
        per_cluster_stats.write(out_line.rstrip())

        # Track number of clusters containing a certain number of KOGs
        try:
            kogs_per_cluster[len(kog_counts)] += 1
        except KeyError:
            kogs_per_cluster[len(kog_counts)] = 1

        # Track number of KOGs spread across a certain number of clusters
        for kog in kog_counts.keys():
            try:
                cluster_counts_by_kog[kog] += 1
            except KeyError:
                cluster_counts_by_kog[kog] = 1

    per_cluster_stats.close()

    # Transform per-KOG cluster counts
    for count in cluster_counts_by_kog.itervalues():
        try:
            clusters_per_kog[count] += 1
        except KeyError:
            clusters_per_kog[count] = 1

    return kogs_per_cluster, clusters_per_kog


#def update_cluster_df(cluster_df, kogs_per_cluster):
#    pass
#
#
#def update_kog_df(kog_df, clusters_per_kog):
#    pass
#
#
#def print_cluster_summary(prefix, kogs_per_cluster):
#    file_name = prefix+"_kogs_per_cluster_summary.Rtab"
#    cluster_summary = open(file_name, 'w')
#
#
#def print_kog_summary(prefix, clusters_per_kog):
#    file_name = prefix+"_clusters_per_kog_summary.Rtab"
#    kog_summary = open(file_name, 'w')


def print_kpc(
        kpc_handle, ordr, frag, ctof, norm, mtrc, infl, kogs_per_cluster):
    """
    """
    for kog_count in sorted(kogs_per_cluster.keys()):
        out_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                   ordr, frag, ctof, norm, mtrc, infl, kog_count,
                   kogs_per_cluster[kog_count])
        kpc_handle.write(out_line)


def print_cpk(
        cpk_handle, ordr, frag, ctof, norm, mtrc, infl, clusters_per_kog):
    """
    """
    for cluster_count in sorted(clusters_per_kog.keys()):
        out_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                   ordr, frag, ctof, norm, mtrc, infl, cluster_count,
                   clusters_per_kog[cluster_count])
        cpk_handle.write(out_line)


if __name__ == "__main__":
    sys.exit(main())
