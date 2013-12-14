#!/usr/bin/env python

"""
This code is unapologetically messy.  It's getting late and this is just a 
small part of what I want to accomplish today.
"""

import sys
from Bio import SeqIO

handles = dict()

for record in SeqIO.parse(open(sys.argv[1], 'r'), 'fasta'):
  org = str(record.id).split('|')[0]

  if not org in handles:
    handles[org] = [open('kog_{0}_{1}.fasta'.format(org,x+1), 'w') for x in range(4)]
#    handles[org] = [open('kog_'+org+'_1.fasta', 'w'),
#                    open('kog_'+org+'_2.fasta', 'w'),
#                    open('kog_'+org+'_3.fasta', 'w')]

  seqID = str(record.id)
  desc = '\t'.join(str(record.description).split('\t')[1:])
  seq = str(record.seq)
  seqLen = len(str(record.seq))

  # Full-length sequences
  handles[org][0].write('>{0}\t{1}\n{2}\n'.format(seqID, desc, seq))

  # Halved sequences
  handles[org][1].write('>{0}\t{1}\n{2}\n'.format(seqID+'---1of2', desc, seq[:int(seqLen/2)]))
  handles[org][1].write('>{0}\t{1}\n{2}\n'.format(seqID+'---2of2', desc, seq[int(seqLen/2):]))

  # Thirded sequences (I hate that work, for the record)
  handles[org][2].write('>{0}\t{1}\n{2}\n'.format(seqID+'---1of3', desc, seq[:int(seqLen/3)]))
  handles[org][2].write('>{0}\t{1}\n{2}\n'.format(seqID+'---2of3', desc, seq[int(seqLen/3):2*int(seqLen/3)]))
  handles[org][2].write('>{0}\t{1}\n{2}\n'.format(seqID+'---3of3', desc, seq[2*int(seqLen/3):]))

  # Quartered sequences
  handles[org][3].write('>{0}\t{1}\n{2}\n'.format(seqID+'---1of4', desc, seq[:int(seqLen/4)]))
  handles[org][3].write('>{0}\t{1}\n{2}\n'.format(seqID+'---2of4', desc, seq[int(seqLen/4):2*int(seqLen/4)]))
  handles[org][3].write('>{0}\t{1}\n{2}\n'.format(seqID+'---3of4', desc, seq[2*int(seqLen/4):3*int(seqLen/4)]))
  handles[org][3].write('>{0}\t{1}\n{2}\n'.format(seqID+'---4of4', desc, seq[3*int(seqLen/4):]))

for org in handles.keys():
  for i in range(3):
    handles[org][i].close()
