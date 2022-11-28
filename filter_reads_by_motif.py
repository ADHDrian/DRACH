import os
import pandas as pd
import pysam
import h5py
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

HOME = os.environ['HOME']

def normalize_signal(input_read):
    s_offset = input_read.attrs["offset"]
    s_range = input_read.attrs["range"]
    s_dig = input_read.attrs["digitisation"]
    signal = (input_read["Dacs"][:] + s_offset) * s_range / s_dig
    shift = input_read.attrs["shift_frompA"]  # median
    scale = input_read.attrs["scale_frompA"]  # mad
    # med = np.median(signal)
    # mad = offset * np.median(abs(signal - med))
    signal = (signal - shift) / scale
    return signal

drach_dir = os.path.join(HOME, 'Data/DRACH')

# base_dir = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293A_wildtype')
base_dir = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293_IVT_2')

taiyaki_dir = os.path.join(base_dir, 'taiyaki')
minimap_dir = os.path.join(base_dir, 'minimap')
outdir = os.path.join(base_dir, 'rodan')
os.makedirs(outdir, exist_ok=True)
motif = 'GGACT'
barcode_len = 17

bed_file = os.path.join(drach_dir, 'motif_{}.bed'.format(motif))
genome_bam_file = os.path.join(minimap_dir, 'aligned_genome.bam.sorted')
taiyaki_file = os.path.join(taiyaki_dir, 'mapped_reads/{}_merge.hdf5'.format(motif))
outfile = os.path.join(outdir, '{}_filtered.hdf5'.format(motif))

num_to_alpha = {
    0 : 'A',
    1 : 'C',
    2 : 'G',
    3 : 'T'
}
alpha_to_num = {v: (k+1) for (k, v) in num_to_alpha.items()}

taiyaki = h5py.File(taiyaki_file, 'r')
bed = pd.read_csv(bed_file, sep='\t', names=['chromosome', 'start_pos', 'end_pos', 'authors', 'score', 'strand', 'motif', 'barcode'])
genome_alignment = pysam.AlignmentFile(genome_bam_file, 'rb')
genome_index = pysam.IndexedReads(genome_alignment)
genome_index.build()

### build read -> barcode dictionary ###
read_ids = list(taiyaki['read_ids'])
dict_id_barcode = {id: [] for id in read_ids}
for qname in tqdm(read_ids):
    entries = genome_index.find(qname)
    for entry in entries:
        chromosome = entry.reference_name
        start_pos = entry.reference_start
        end_pos = entry.reference_end
        bed_chromosome = bed[bed['chromosome']==chromosome]
        matches = bed_chromosome[(bed_chromosome['start_pos']>=start_pos) & (bed_chromosome['end_pos']<=end_pos)]
        if len(matches)==0:
            # print('No matches found for {}'.format(qname))
            continue
        for _, row in matches.iterrows():
            dict_id_barcode[qname].append(row['barcode'])

### cut out barcoded-sections in each read ###
chunk_len = 1024
min_spb = 20
max_spb = 70

id_signal_label = []
for read_id in tqdm(read_ids):
    barcodes = dict_id_barcode[read_id]
    if len(barcodes)==0:
        # print('No barcodes found for {}'.format(read_id))
        continue

    this_read = taiyaki['Reads'][read_id]
    # dacs = np.array(this_read['Dacs'])
    signal = normalize_signal(this_read)
    ref_to_signal = np.array(this_read['Ref_to_signal'])
    reference = np.array(this_read['Reference'])

    label = ''.join([num_to_alpha[i] for i in reference])

    for bc in barcodes:
        bc_start = label.find(bc[::-1])
        if (bc_start < 0):
            # print('Barcode not found. Skipping...')
            continue

        signal_start = ref_to_signal[bc_start]
        signal_end = signal_start + chunk_len
        if signal_end>=len(signal):
            continue
        bc_signal = signal[signal_start:signal_end]
        bc_end = np.searchsorted(ref_to_signal, signal_end)
        if bc_end>=len(label):
            continue
        bc_reference = reference[bc_start:bc_end]

        num_bases = len(bc_reference)
        spb = chunk_len / num_bases

        if (len(bc_signal)<chunk_len) or (num_bases<(barcode_len//2+1)) or (spb<min_spb) or (spb>max_spb):
            # print('Num. bases not matching requirements')
            continue

        # bc_end = bc_start + len(bc)
        # bc_reference = reference[bc_start:bc_end]
        # bc_ref_to_signal = ref_to_signal[bc_start:bc_end+1]
        # bc_base_loc = (bc_ref_to_signal[1:] + bc_ref_to_signal[:-1]) / 2

        # if (bc_start<0) or ((bc_ref_to_signal[-1]-bc_ref_to_signal[0])>chunk_len):
        #     # print('Skipping {}...'.format(read_id))
        #     continue

        id_signal_label.append((read_id, bc_signal, bc_reference))

    ### debug ###
    # plt.figure(figsize=(10, 5))
    # plt.plot(signal)
    # for loc in motif_ref_to_signal:
    #     plt.axvline(loc, c='r')
    # plt.xlim([motif_ref_to_signal[0], motif_ref_to_signal[-1]])
    # # plt.xticks(motif_base_loc, motif[::-1])
    # plt.xticks([motif_base_loc[10]], [motif[10]], c='r')
    # plt.title(read_id)
    # plt.savefig(os.path.join(outdir, 'motif_selection_{}.png'.format(read_id)), bbox_inches='tight')
    # plt.close('all')

### filter spb ###
# spbs = [len(sig)/len(la) for (_, sig, la) in id_signal_label]
# filtered_id_signal_label = []
# for id_sig_lab, spb in zip(id_signal_label, spbs):
#     if (spb>=min_spb) and (spb<=max_spb):
#         filtered_id_signal_label.append(id_sig_lab)

### write out events ###
num_events = len(id_signal_label)
labels_len = [len(label) for (_, signal, label) in id_signal_label]
max_label_len = max(labels_len)
# rev_motif_num = np.array([alpha_to_num[x] for x in motif[::-1]])
# label_len = len(rev_motif_num)
# id_len = len(id_signal_label[0][0])

with h5py.File(outfile, "w") as h5_out:
    h5_out_read_id = h5_out.create_dataset("read_ids", dtype=h5py.string_dtype(encoding='ascii'), shape=(num_events,))
    h5_out_event = h5_out.create_dataset("events", dtype="float32", shape=(num_events, chunk_len))
    h5_out_label = h5_out.create_dataset("labels", dtype="int64", shape=(num_events, max_label_len))
    h5_out_label_len = h5_out.create_dataset("labels_len", dtype="int64", shape=(num_events,))

    for i, (read_id, signal, label) in enumerate(id_signal_label):
        # pad_left = (chunk_len - len(signal)) // 2
        # h5_out_event[i] = np.pad(signal, (pad_left, chunk_len - pad_left - len(signal)))
        h5_out_read_id[i] = read_id
        h5_out_event[i] = signal
        h5_out_label[i] = np.pad(label+1, (0, max_label_len - len(label)))
        h5_out_label_len[i] = labels_len[i]