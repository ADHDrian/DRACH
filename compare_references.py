import pandas as pd
from Bio import SeqIO, Seq
from collections import Counter

base_ref_file = '/home/adrian/Data/transcriptomes/human_reference.fasta'
compare_ref_file = '/home/adrian/Data/transcriptomes/gencode.v42.transcripts.fa'

base_tx = {}
with open(base_ref_file) as handle:
    for record in SeqIO.parse(handle, format='fasta'):
        base_tx[record.id.split('.')[0]] = str(record.seq)

compare_tx = {}
with open(compare_ref_file) as handle:
    for record in SeqIO.parse(handle, format='fasta'):
        compare_tx[record.id.split('.')[0]] = str(record.seq)

common_keys = set(base_tx.keys()).intersection(set(compare_tx.keys()))
print('{} in base, {} in compare'.format(len(base_tx), len(compare_tx)))
print('{} in common'.format(len(common_keys)))

missing_keys = set(base_tx.keys()) - set(compare_tx.keys())
print('{} missing from base reference'.format(len(missing_keys)))

missing_tags = []
with open(base_ref_file) as handle:
    for record in SeqIO.parse(handle, format='fasta'):
        id = record.id
        if id.split('.')[0] in missing_keys:
            missing_tags.append(id.split('|')[-2])

for k, v in Counter(missing_tags).most_common():
    print('{}: {}'.format(k, v))