import os
import pandas as pd
import pysam
from Bio import SeqIO, Align
# from Bio.Align import substitution_matrices
# aligner = Align.PairwiseAligner()
# aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

from collections import Counter

# BARCODE_SPAN = 8
# WORKSPACE = '/prj/Isabel_IVT_Nanopore/HEK293A_wildtype/Jessica_HEK293/HEK293A_2/20190409_1503_GA10000_FAK10978_2e75d7be/achan/taiyaki'
# MOTIF = 'GGACT'
# miclip_file = '/home/achan/Data/DRACH/miCLIP_union_flat_exclude_Y_chromosome.bed'
# genome_ref = '/home/achan/Data/genomes/GRCh38_96.fa'
# genome_bam_file = os.path.join(WORKSPACE, 'aligned_genome.bam.sorted')
# transcriptome_bam_file = os.path.join(WORKSPACE, 'aligned_transcriptome.bam.sorted')
# out_transcriptome_bam_file = os.path.join(WORKSPACE, 'transcriptome_{}.bam'.format(MOTIF))

def parse_genome_ref(genome_ref):
    genome = {}
    with open(genome_ref, 'r') as handle:
        for record in SeqIO.parse(handle, format='fasta'):
            id = record.id
            if id.isnumeric() or id in ['X', 'Y']:
                print(id)
                genome[id] = record.seq
    return genome

def parse_miclip_file(miclip_file, genome, motif):
    miclip = pd.read_csv(miclip_file, sep='\t', names=['chromosome', 'start_pos', 'end_pos', 'authors', 'score', 'strand'])
    mask_chrom = (pd.to_numeric(miclip['chromosome'], errors='coerce').notnull() | miclip['chromosome'].isin(['X', 'Y']))
    miclip = miclip[mask_chrom]
    miclip = miclip.reset_index(drop=True)

    barcodes = []
    motifs = []
    for ind, row in miclip.iterrows():
        chromosome = str(row['chromosome'])
        strand = row['strand']
        pos = row['start_pos']
        barcode = genome[chromosome][pos-BARCODE_SPAN:pos+BARCODE_SPAN+1] if strand=='+' else genome[chromosome][pos-BARCODE_SPAN:pos+BARCODE_SPAN+1].reverse_complement()
        motif = barcode[BARCODE_SPAN-2:BARCODE_SPAN+3]
        # print(ind, strand, motif)
        barcodes.append(str(barcode))
        motifs.append(str(motif))
    miclip['barcode'] = barcodes
    miclip['motif'] = motifs

    most_common_motifs = Counter(motifs).most_common(50)
    for k, v in most_common_motifs:
        print(k, v)

    df_motif = miclip[miclip['motif'] == motif]
    df_motif.reset_index(drop=True, inplace=True)
    df_motif.to_csv(os.path.join(WORKSPACE, 'motif_{}.bed'.format(motif)), sep='\t', header=None, index=None)

    return df_motif

def parse_genome_alignment(genome_bam_file, df_motif):
    genome_bam = pysam.AlignmentFile(genome_bam_file, 'rb')
    site_ind_read_ids = {}
    for site_ind, this_site in df_motif.iterrows():
        # print(this_site)
        this_site_fetch = genome_bam.fetch(this_site.chromosome, this_site.start_pos, this_site.end_pos)
        read_ids = []
        for this_read in this_site_fetch:
            # print(str(this_read).split('\t'))
            qname = this_read.qname
            read_ids.append(qname)
        site_ind_read_ids[site_ind] = read_ids
    unique_read_ids = list(set([x for l in site_ind_read_ids.values() for x in l]))
    return unique_read_ids

def write_motif_filtered_transcriptome_alignment(transcriptome_bam_file, unique_read_ids):
    transcriptome_bam = pysam.AlignmentFile(transcriptome_bam_file, 'rb')
    transcriptome_indexed = pysam.IndexedReads(transcriptome_bam)
    transcriptome_indexed.build()
    header = transcriptome_bam.header.copy()
    out = pysam.AlignmentFile(out_transcriptome_bam_file, 'wb', header=header)
    for this_read_id in unique_read_ids:
        try:
            transcriptome_indexed.find(this_read_id)
        except KeyError:
            pass
        else:
            iterator = transcriptome_indexed.find(this_read_id)
            for x in iterator:
                out.write(x)
    out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Data.')
    parser.add_argument("--workspace", default=None, type=str)
    parser.add_argument("--motif", default=None, type=str)
    parser.add_argument("--span", default=8, type=int, help="Architecture file")
    parser.add_argument("--ref", default=None, type=str)
    parser.add_argument("--miclip", default=None, type=str)
    parser.add_argument("--genome_bam", default=None, type=str)
    parser.add_argument("--transcriptome_bam", default=None, type=str)
    parser.add_argument("--outfile", default=None, type=str)
    args = parser.parse_args()

    BARCODE_SPAN = args.span

    genome = parse_genome_ref(args.ref)
    df_motif = parse_miclip_file(args.miclip, genome, args.motif)
    read_ids = parse_genome_alignment(args.genome_bam, df_motif)
    write_motif_filtered_transcriptome_alignment(args.transcriptome_bam, args.outfile, read_ids)