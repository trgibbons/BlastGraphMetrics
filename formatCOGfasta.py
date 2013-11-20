#!/usr/bin/env python

import sys
from Bio import SeqIO


whog = open(sys.argv[1], 'r')
fasta = open(sys.argv[2], 'r')
out = open(sys.argv[3], 'w')

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


for record in SeqIO.parse(fasta, 'fasta'):
  old_id = record.id
  seq = record.seq
  if old_id in new_ids:
    out.write('>{0}\n{1}\n'.format(new_ids[old_id], seq))


whog.close()
fasta.close()
out.close()
