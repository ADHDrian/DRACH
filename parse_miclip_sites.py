import os
import pandas as pd
import pysam
from Bio import SeqIO, Align
# from Bio.Align import substitution_matrices
# aligner = Align.PairwiseAligner()
# aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

from collections import Counter

miclip_file = '/home/achan/Data/DRACH/miCLIP_union_flat_exclude_Y_chromosome.bed'
genome_file = '/home/achan/Data/genomes/GRCh38_96.fa'
outdir = '/home/achan/Data/DRACH'

miclip = pd.read_csv(miclip_file, sep='\t', names=['chromosome', 'start_pos', 'end_pos', 'authors', 'score', 'strand'])
mask_chrom = (pd.to_numeric(miclip['chromosome'], errors='coerce').notnull() | miclip['chromosome'].isin(['X', 'Y']))
miclip = miclip[mask_chrom]
miclip = miclip.reset_index(drop=True)

genome = {}
with open(genome_file, 'r') as handle:
    for record in SeqIO.parse(handle, format='fasta'):
        id = record.id
        if id.isnumeric() or id in ['X', 'Y']:
            print(id)
            genome[id] = record.seq

motifs = []
SPAN = 10
for ind, row in miclip.iterrows():
    chromosome = str(row['chromosome'])
    strand = row['strand']
    pos = row['start_pos']
    motif = genome[chromosome][pos-SPAN:pos+SPAN+1] if strand=='+' else genome[chromosome][pos-SPAN:pos+SPAN+1].reverse_complement()
    print(ind, strand, motif)
    motifs.append(str(motif))

miclip['motif'] = motifs

most_common_motifs = Counter(motifs).most_common(50)
for k, v in most_common_motifs:
    print(k, v)

### choose motif ###
this_motif = most_common_motifs[0][0]
df_motif = miclip[miclip['motif'] == this_motif]
df_motif['start_pos'] = df_motif['start_pos'] - SPAN
df_motif['end_pos'] = df_motif['end_pos'] + SPAN
df_motif.loc[:, df_motif.columns!='motif'].to_csv(os.path.join(outdir, 'motif_{}.bed'.format(this_motif)), sep='\t', header=None, index=None)

### parse filtered bam file ###
workspace = '/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/achan/taiyaki'
genome_bam_file = os.path.join(workspace, 'aligned_genome.bam.sorted')
genome_bam = pysam.AlignmentFile(genome_bam_file, 'rb')

transcriptome_bam_file = os.path.join(workspace, 'aligned_transcriptome.bam.sorted')
transcriptome_bam = pysam.AlignmentFile(transcriptome_bam_file, 'rb')

id_seq = []
for _, this_site in df_motif.iterrows():
    # print(this_site)
    this_site_iter = bam.fetch(this_site.chromosome, this_site.start_pos, this_site.end_pos)
    for this_read in this_site_iter:
        # print(str(this_read).split('\t'))
        qname = this_read.qname
        pos = this_read.pos
        seq = this_read.seq
        site_read_ids.append(qname)


rev_motif = this_motif[::-1]
with open(os.path.join(workspace, 'read_references_{}.fasta'.format(this_motif)), 'w') as outfile:
    for this_read_id in site_read_ids:
        outfile.write('>{}\n'.format(this_read_id))
        outfile.write('{}\n'.format(rev_motif))