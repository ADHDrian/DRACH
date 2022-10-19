import pandas as pd
from Bio import SeqIO

transcript_file = '/home/adrian/Data/transcriptomes/Homo_sapiens.GRCh38.cdna.all.fa'
mod_file = '/home/adrian/Data/Nanopore_HEK293/miCLIP2/miCLIP_m6A_sites_filtered_transcriptome.tsv'
motif_file = '/home/adrian/Data/Nanopore_HEK293/miCLIP2/drach_motifs.fa'

mod_table = pd.read_csv(mod_file, sep='\t', usecols=['group', 'end', 'tx_id'])

transcript = {}
with open(transcript_file) as handle:
    for record in SeqIO.parse(handle, format='fasta'):
        transcript[record.id.split('.')[0]] = str(record.seq)

SPAN = 10
drach_motifs = {}
for ind, row in mod_table.iterrows():
    if row.tx_id not in transcript.keys():
        continue
    seq = transcript[row.tx_id]
    mod_pos = row.end - 1
    if (mod_pos>(len(seq)-SPAN)) or (mod_pos<SPAN):
        continue

    motif = seq[mod_pos-SPAN:mod_pos+SPAN+1]
    drach_motifs[row.tx_id] = motif

    if (
            (motif[10-2] in ['A', 'G', 'T']) &
            (motif[10-1] in ['A', 'G']) &
            (motif[10]=='A') &
            (motif[10+1]=='C') &
            (motif[10+2] in ['A', 'C', 'T'])
    ): print('DRACH!')
    else: print('NOT!')