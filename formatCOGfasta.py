#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO


def get_parsed_args():
  """Parse the command line arguments

  Parses command line arguments using the argparse package, which is a standard
  Python module starting with version 2.7.

  Args:
    None, argparse fetches them from user input

  Returns:
    args: An argparse.Namespace object containing the parsed arguments
  """
  parser = argparse.ArgumentParser(
             description='Generate a set of graphs from a tab-delimited '+
                         'BLASTP or BLASTN file such that the first two '+
                         'columns contain the query and subject IDs, '+
                         'respectively, and the last four columns contain, in '+
                         'order: E-value, bit score, query length, subject '+
                         'length')

  parser.add_argument('--whog', dest='whog', default=False,
                      type=argparse.FileType('r'),
                      help='whog file that comes with the COG/COG database')
  parser.add_argument('--twog', dest='twog', default=False,
                      type=argparse.FileType('r'),
                      help='twog file that comes with the COG/KOG database')
  parser.add_argument('fasta', type=argparse.FileType('r'),
                      help='Either the myva (COG) or kyva (KOG) FASTA file')
  parser.add_argument('out', type=argparse.FileType('w'),
                      help='Output file name')

  args = parser.parse_args()

  return args



def parse_whog(whog):
  """Prepare new COG/COG seq IDs from the whog file"""
  new_ids = dict()

  for line in whog:
    temp = line.strip().split()
  
    if line == '\n':
      continue
    elif line == '_______\n':
      continue
  
    elif temp[0][0] == '[':
      fun = str(temp[0])
      cog = str(temp[1])
      des = ' '.join(temp[2:])
  
    else:
      org = temp[0].rstrip(':')
      for seq in temp[1:]:
        new_ids[seq] = '{0}|{1}___{2}\t{3}\t{4}'.format(org, seq, cog, fun, des)

  return new_ids



def parse_twog(twog):
  """Prepare new KOG seq IDs from the twog file"""




def reformat_fasta(new_ids, fasta, out):
  for record in SeqIO.parse(fasta, 'fasta'):
    old_id = record.id
    seq = record.seq
    if old_id in new_ids:
      out.write('>{0}\n{1}\n'.format(new_ids[old_id], seq))



def main(argv=None):
  """Where the magic happens!
  
  The main() function coordinates calls to all of the other functions in this 
  program in the hope that, by their powers combined, useful work will be done.
  
  Args:
    None

  Returns:
    An exit status (hopefully 0)

  Raises:
    SystemExit: Thrown when incompatible options are selected
  """
  if argv == None:
    argv = sys.argv

  args = get_parsed_args()

  if args.whog and args.twog:
    raise SystemExit('The --whog and --twog arguments are incompatible. Please '+
                    'specify only one and provide the corresponding FASTA '+
                    'file (myva for whog, kyva for twog)')
  elif args.whog:
    new_ids = parse_whog(args.whog)
  elif args.twog:
    new_ids = parse_twog(args.twog)
  else:
    raise SystemExit('Must provide exactly one of --whog and --twog along with '+
                    'the corresponding FASTA file (myva for whog, kyva for '+
                    'twog)')

  reformat_fasta(new_ids=new_ids, fasta=args.fasta, out=args.out)



if __name__ == "__main__":
  sys.exit(main())

