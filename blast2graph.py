#!/usr/bin/env python

# Standard Python libraries
import sys
import os
import argparse

# Third-party libraries
import networkx as nx



def catch_argument_errors(args):
  """Attempt to catch argument and file format errors

  Args:
    args: The argparse.Namespace object returned by get_parse_args

  Returns:
    None (hopefully)

  Raises:
    TODO: Implement proper Python exception handling
  """
  kill_switch = 0

  if os.fstat(args.both.fileno()).st_size == 0:
    args.both = False
    if not len(args.self)>0 and len(args.cross)>0:
      sys.stderr.write(
        'Must have at least one file flagged as "both", or at least one each '+
        'flagged as "self" and "cross", repectively.\n')
      kill_switch += 1

  if kill_switch > 0:
    sys.exit()



def get_parsed_args():
  """Parse the command line arguments

  Parses command line arguments using the argparse package, which is a standard
  Python module starting with version 2.7.

  Args:
    None, argparse fetches them from user input

  Returns:
    args: An argparse.Namespace object containing the parsed arguments

  Raises:
    None
  """
  parser = argparse.ArgumentParser(
             description='Compute a set of Bit Score Ratios from either '+
                         'BLASTP or BLASTN hits')

  # Group: IO options
  parser.add_argument('-s', '--self', dest='self', nargs='+',
                      type=argparse.FileType('r'),
                      help='Tab-delimited BLAST file(s) containing self '+
                           'alignments (non-self alignments will be ignored)')
  parser.add_argument('-c', '--cross', dest='cross', nargs='+',
                      type=argparse.FileType('r'),
                      help='Tab-delimited BLAST file(s) containing non-self '+
                           'alignments (self alignments will be ignored)')
  parser.add_argument('-b', '--both', dest='both', nargs='?',
                      type=argparse.FileType('r'), default=sys.stdin,
                      help='Tab-delimited BLAST file containing a combination '+
                           'of self and non-self alignments [def=stdin]')
  parser.add_argument('-o', '--out', dest='out', nargs='?',
                      type=argparse.FileType('w'), default=sys.stdout,
                      help='Name for MCL-formatted graph output file '+
                           '[def=stdout]')
  parser.add_argument('-d', '--debug', dest='debug',
                      type=argparse.FileType('w'), default=False,
                      help='File for capturing debug messages')

  # Group: Formatting options
  parser.add_argument('--bscol', dest='bscol',
                      action='store', type=int, default=12,
                      help='One-indexed column containing pairwise bit scores '+
                           '(not required if files include standard header '+
                           'lines) [def=12]')
  parser.add_argument('--idchar', dest='idchar', action='store', default='|',
                      help='The character used to separate the organism ID '+
                           'from the rest of the sequence header [def="|"]')
  parser.add_argument('--idlist', dest='idlist', action='store', default=None,
                      help='Text file containing a list of organisms IDs, '+
                           'which must occur at the beginning of each '+
                           'sequence ID ("--idchar" option must also be '+
                           'provided if any IDs are prefixes of other IDs)')

  # Group: Behavioral options
  parser.add_argument('-n', '--normalize', dest='norm',
                      action='store_true', default=False,
                      help='Normalize edge weights according to intra- and '+
                           'inter-organism averages (requires the "-i" option)')
  parser.add_argument('-r', '--reciprocal', dest='recip',
                      action='store_true', default=False,
                      help='Divide alignment bit score by max(self1,self2) '+
                      'instead of min(self1,self2) for increased stringency '+
                      '[def=False]')

  # Group: TODO
  parser.add_argument('-t', '--translated', action='store_true', default=False,
                      help='One or more of the data sets were '+
                           'bioinformatically translated by BLAST (requires '+
                           'the "--idchar" and/or "--idlist" option(s) be '+
                           'used with a file containing a second column '+
                           'specifying either "nuc" or "pro" for each '+
                           'organism ID listed in the first column')
  parser.add_argument('-m', '--merge', dest='merge',
                      action='store_true', default=False,
                      help='Merge sequences from a single organism when they '+
                           'have non-overlapping alignments to the same '+
                           'target sequence (helpful for highly-fragmented '+
                           'assemblies, requires the "-i" option)')

  args = parser.parse_args()

  catch_argument_errors(args)

  return args



def get_self_bit_scores(
      bsr_graph, self_handle, bscol=11, norm=False, idchar=None, org_ids=None):
  """Get bit scores from full-length self-alignments

  Searches an open file for tab-delimited BLAST hit records where the query and
  reference IDs are identical.  It is important that BLAST is run with option
  "--soft_masking true" or the self-alignments are unlikely to be full length.

  BLAST does not guarantee particular boundaries for a given alignment, and 
  reciprocal alignments frequently have slightly different boundries and 
  scores. Neither is more or less valid, so I see no reason to not use the 
  larger of the two. I therefore check to make sure subsequent hits do not have 
  a greater score than the ones that have already been found. One consequence 
  of this is that the intra- and inter-organism average scores can not be 
  calculate on the fly.

  Args:
    bsr_graph: A NetworkX graph data structure (does not need to be empty)
    self_handle: An open file handle containing self-alignments (can contain
      other alignments and/or comment lines beginning with a hash '#' character)
    bscol: Column containing the bit scores

  Returns:
    Nothing, all data structures are edited in place

  Raises:
    None
  """
  for line in self_handle:
    temp = line.strip().split()
    if temp[0][0] == "#":
      continue

    if temp[0] == temp[1]:
      seq_id = str(temp[0])
      bitscore = float(temp[bscol])

      if norm:
        org_ids.add(seq_id.split(idchar)[0])

      if not bsr_graph.has_node(seq_id):
        bsr_graph.add_node(seq_id, sbs=bitscore, org=None)
      elif bitscore > bsr_graph.node[seq_id]['sbs']:
        bsr_graph.node[seq_id]['sbs'] = bitscore



def get_cross_bit_scores(bsr_graph, cross_handle, bscol=11, recip=False):
  """Get bit scores from full-length alignments between different sequences

  Searches an open file for tab-delimited BLAST hit records where the query and
  reference IDs are not identical. It is important that BLAST is run with option
  "--soft_masking true" or the self-alignments are unlikely to be full length.

  BLAST does not guarantee particular boundaries for a given alignment, and 
  reciprocal alignments frequently have slightly different boundries and 
  scores. Neither is more or less valid, so I see no reason to not use the 
  larger of the two. I therefore check to make sure subsequent hits do not have 
  a greater score than the ones that have already been found. One consequence 
  of this is that the intra- and inter-organism average scores can not be 
  calculate on the fly.

  Args:
    bsr_graph: A NetworkX graph data structure containing self-alignment scores
    self_handle: An open file handle containing non-self-alignments (can contain
      other alignments and/or comment lines beginning with a hash '#' character)
    bscol: Column containing the bit scores

  Returns:
    Nothing, all data structures are edited in place

  Raises:
    None
  """
  for line in cross_handle:
    temp = line.strip().split()
    if temp[0][0] == "#":
      continue

    if temp[0] != temp[1]:
      qry_id = str(temp[0])
      ref_id = str(temp[1])
      bitscore = float(temp[bscol])

      if bsr_graph.has_node(qry_id) and bsr_graph.has_node(ref_id):
        qry_sbs = bsr_graph.node[qry_id]['sbs']
        ref_sbs = bsr_graph.node[ref_id]['sbs']

        if recip:
          bsr = bitscore / max(qry_sbs, ref_sbs)
        else:
          bsr = bitscore / min(qry_sbs, ref_sbs)

        if not bsr_graph.has_edge(qry_id, ref_id):
          bsr_graph.add_edge(qry_id, ref_id, bsr=bsr)
        elif bsr > bsr_graph[qry_id][ref_id]['bsr']:
          bsr_graph[qry_id][ref_id]['bsr'] = bsr



def compute_organism_averages(bsr_graph, org_ids, idchar):
  """Compute average scores between and within each pair of organisms
  
  Args:
    bsr_graph: A NetworkX graph data structure containing self-alignment scores
    org_ids: A set containing each organism ID
  Returns:
    org_avgs: A NetworkX graph data structure containing the total number of 
      edges between each pair of organisms, the combined
  """
  org_avgs = nx.Graph()

  for qry_id, ref_id, edata in bsr_graph.edges(data=True):
    qry_org = qry_id.split(idchar)[0]
    ref_org = ref_id.split(idchar)[0]
    if org_avgs.has_edge(qry_org, ref_org):
      org_avgs[qry_org][ref_org]['cnt'] += 1
      org_avgs[qry_org][ref_org]['sum'] += edata['bsr']
    else:
      org_avgs.add_edge(qry_org, ref_org, cnt=1, sum=edata['bsr'], avg=None)

  org_avgs.add_node('global', cnt=0, sum=0.)
  for qry_org, ref_org, edata in org_avgs.edges(data=True):
    org_avgs.node['global']['cnt'] += edata['cnt']
    org_avgs.node['global']['sum'] += edata['sum'] #float
    org_avgs[qry_org][ref_org]['avg'] = edata['sum']/edata['cnt']

  glb_cnt = org_avgs.node['global']['cnt']
  glb_sum = org_avgs.node['global']['sum'] #float
  org_avgs.node['global']['avg'] = glb_sum/glb_cnt

  return org_avgs



def normalize_bit_score_ratios(bsr_graph, org_avgs, idchar):
  """Convert Bit Scores into Bit Score Ratios

  Iterates through the edges in a NetworkX graph, dividing all cross-alignment
  scores by either the smaller or larger of the two self-alignment scores

  Convert Bit Scores into Bit Score Ratios and account for intra-/inter- 
  organism differences, if requested
  """
  glb_avg = org_avgs.node['global']['avg']
  for qry_id, ref_id, edata in bsr_graph.edges(data=True):
    qry_org = qry_id.split(idchar)[0]
    ref_org = ref_id.split(idchar)[0]
    scale = glb_avg / org_avgs[qry_org][ref_org]['avg']
    raw_avg = bsr_graph[qry_id][ref_id]['bsr']
    bsr_graph[qry_id][ref_id]['bsr'] = raw_avg * scale



def print_mcl_input_file(bsr_graph, out_handle):
  """Print an MCL-formatted graph file"""
  for line in nx.generate_edgelist(bsr_graph, delimiter='\t', data=['bsr']):
    out_handle.write(line+'\n')



def debug_print(debug_handle, org_avgs=None, bsr_graph=None):
  """Dump some data structures for inspection"""
  if org_avgs:
    for u, v, data in bsr_graph.edges(data=True):
      debug_handle.write(str(u))
  for node, data in bsr_graph.nodes(data=True):
    debug_handle.write(str(node)+"\t"+repr(data)+"\n")



def main(argv=None):
  """Where the magic happens!
  
  The main() function coordinates calls to all of the other functions in this 
  program in the hope that, by their powers combined, useful work will be done.
  
  Args:
    None

  Returns:
    An exit status (hopefully 0)

  Raises:
    Hopefully nothing
  """
  if argv == None:
    argv = sys.argv

  args = get_parsed_args() 

  bsr_graph = nx.Graph()
  org_ids = set()


  if args.self:
    for self_handle in args.self:
      get_self_bit_scores(
        bsr_graph=bsr_graph, self_handle=self_handle, bscol=args.bscol-1,
        norm=args.norm, idchar=args.idchar, org_ids=org_ids)

  if args.both:
    get_self_bit_scores(
      bsr_graph=bsr_graph, self_handle=args.both, bscol=args.bscol-1,
      norm=args.norm, idchar=args.idchar, org_ids=org_ids)
    args.both.seek(0)
    get_cross_bit_scores(
      bsr_graph=bsr_graph, cross_handle=args.both, bscol=args.bscol-1,
      recip=args.recip)

  if args.cross:
    for cross_handle in args.cross:
      get_cross_bit_scores(bsr_graph, cross_handle, args.bscol-1)

  if args.norm:
    org_avgs = compute_organism_averages(
                 bsr_graph=bsr_graph, org_ids=org_ids, idchar=args.idchar)
    normalize_bit_score_ratios(
      bsr_graph=bsr_graph, org_avgs=org_avgs, idchar=args.idchar)

  print_mcl_input_file(bsr_graph, args.out)


if __name__ == "__main__":
  sys.exit(main())

