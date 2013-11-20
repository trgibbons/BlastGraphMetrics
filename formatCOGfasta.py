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
                      help='whog file that comes with the COG database')
  parser.add_argument('--kog', dest='kog', default=False,
                      type=argparse.FileType('r'),
                      help='kog file that comes with the KOG database')
  parser.add_argument('fasta', type=argparse.FileType('r'),
                      help='Either the myva (COG) or kyva (KOG) FASTA file')
  parser.add_argument('out', type=argparse.FileType('w'),
                      help='Output file name')

  args = parser.parse_args()

  return args



def parse_whog(whog):
  """Prepare new COG/COG seq IDs from the whog file
  
  Args:
    whog: A readable handle for the whog file
  
  Returns:
    new_ids: A dictionary mapping expanded sequence IDs onto the old ones
  """
  new_ids = dict()

  for line in whog:
    temp = line.strip().split()
  
    if not temp:
      continue

    elif temp[0] == '_______':
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



def parse_kog(kog):
  """Prepare new KOG seq IDs from the kog file
  
  Args:
    kog: A readable handle for the kog file
  
  Returns:
    new_ids: A dictionary mapping expanded sequence IDs onto the old ones
  """
  new_ids = dict()
  debug = open('debug', 'w')

  for line in kog:
    debug.write(line)
    temp = line.strip().split()
    debug.write(','.join(temp)+'\n')
  
    if not temp:
      continue

    elif temp[0][0] == '[':
      fun = str(temp[0])
      cog = str(temp[1])
      des = ' '.join(temp[2:])

    else:
      org = temp[0].rstrip(':')
      for seq in temp[1:]:
        new_ids[seq] = '{0}|{1}___{2}\t{3}\t{4}'.format(org, seq, cog, fun, des)

  debug.close()

  return new_ids



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

  if args.whog and args.kog:
    raise SystemExit('The --kog and --twog arguments are incompatible. Please '+
                    'specify only one and provide the corresponding FASTA '+
                    'file (myva for whog, kyva for kog)')
  elif args.whog:
    new_ids = parse_whog(args.whog)
  elif args.kog:
    new_ids = parse_kog(args.kog)
  else:
    raise SystemExit('Must provide exactly one of --whog and --kog along with '+
                    'the corresponding FASTA file (myva for whog, kyva for '+
                    'kog)')

  reformat_fasta(new_ids=new_ids, fasta=args.fasta, out=args.out)



if __name__ == "__main__":
  sys.exit(main())

