import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

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

taiyaki_file = '/home/adrian/Data/Isabel_IVT_ONT/HEK293A_wildtype/taiyaki/mapped_reads/TAATCCCAGCACTTTGGGAGG_merge.hdf5'
motif = 'TAATCCCAGCACTTTGGGAGG'

outdir = '/home/adrian/Data/Isabel_IVT_ONT/HEK293A_wildtype/taiyaki/mapped_reads'
outfile = os.path.join(outdir, 'TAATCCCAGCACTTTGGGAGG_filtered.hdf5')

num_to_alpha = {
    0 : 'A',
    1 : 'C',
    2 : 'G',
    3 : 'T'
}

alpha_to_num = {v: (k+1) for (k, v) in num_to_alpha.items()}

h5_in = h5py.File(taiyaki_file, 'r')

max_event_len = 1024

id_signals = {}
for read_id in h5_in['read_ids']:
    this_read = h5_in['Reads'][read_id]
    # dacs = np.array(this_read['Dacs'])
    signal = normalize_signal(this_read)
    ref_to_signal = np.array(this_read['Ref_to_signal'])
    reference = np.array(this_read['Reference'])

    label = ''.join([num_to_alpha[i] for i in reference])

    motif_start = label.find(motif[::-1])
    motif_end = motif_start + len(motif)

    motif_reference = reference[motif_start:motif_end]
    motif_ref_to_signal = ref_to_signal[motif_start:motif_end+1]
    motif_base_loc = (motif_ref_to_signal[1:] + motif_ref_to_signal[:-1]) / 2

    if (motif_start<0) or ((motif_ref_to_signal[-1]-motif_ref_to_signal[0])>max_event_len):
        print('Skipping {}...'.format(read_id))
        continue

    id_signals[read_id] = signal[motif_ref_to_signal[0]:motif_ref_to_signal[-1]]

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

num_events = len(id_signals)
rev_motif_num = np.array([alpha_to_num[x] for x in motif[::-1]])
label_len = len(rev_motif_num)

with h5py.File(outfile, "w") as h5_out:
    h5_out_event = h5_out.create_dataset("events", dtype="float32", shape=(num_events, max_event_len))
    h5_out_label = h5_out.create_dataset("labels", dtype="int64", shape=(num_events, len(motif)))
    h5_out_label_len = h5_out.create_dataset("labels_len", dtype="int64", shape=(num_events,))

    for i, signal in enumerate(id_signals.values()):
        h5_out_event[i] = np.pad(signal, (0, 1024 - len(signal)))
        h5_out_label[i] = rev_motif_num
        h5_out_label_len[i] = label_len