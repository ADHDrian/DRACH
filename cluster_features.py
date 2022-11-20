import pickle
from Bio.pairwise2 import align, format_alignment
import numpy as np
import matplotlib.pyplot as plt

motif = 'TAATCCCAGCACTTTGGGAGG'
span = 10
ref = motif[::-1]

def calc_mod_pos(in_str, burner=ref[:(span+1)]):
    pos = 0
    while len(burner)>0:
        if in_str[0]=='-':
            pos += 1
            in_str = in_str[1:]
        elif in_str[0]==burner[0]:
            pos += 1
            in_str = in_str[1:]
            burner = burner[1:]
        else:
            print('WTH is this?', in_str[0])
    return pos - 1


file = '/home/adrian/Data/DRACH/TAATCCCAGCACTTTGGGAGG/loss_pred_label_feature.pkl'

vocab = { 1:"A", 2:"C", 3:"G", 4:"T" }

with open(file, 'rb') as f:
    loss_pred_label_feature = pickle.load(f)

mat_feature = []
for event in loss_pred_label_feature.values():
    loss = event['loss']
    pred_label = event['pred_label']
    feature = event['feature']

    alignments = align.globalxx(ref, pred_label)
    aligned = alignments[0]

    if aligned.score >= 15:
        ref_mod_pos = calc_mod_pos(aligned.seqA)
        pred_mod_pos = ref_mod_pos - len([x for x in aligned.seqB[:ref_mod_pos] if x not in vocab.values()])

        # print(format_alignment(*aligned))
        # print(aligned.seqB[pred_mod_pos-1:pred_mod_pos+2])
        print(aligned.seqA[:ref_mod_pos+1])
        print(pred_label[:pred_mod_pos+1])
        print('\n')

        mod_feature = feature[pred_mod_pos]
        mat_feature.append(mod_feature)

        # plt.figure()
        # plt.plot(mod_feature)
        # plt.show()
    else:
        print('Score {} too low'.format(aligned.score))
        print('\n')
mat_feature = np.vstack(mat_feature)

### calculate pairwise distance ###
from scipy.spatial.distance import pdist
c_vec = 1.0 - pdist(mat_feature, 'cosine')

plt.figure(figsize=(10, 6))
plt.hist(c_vec, bins=50, range=[-1, 1])
plt.xlabel('Cosine similarity')
plt.ylabel('Counts')
plt.title('{}\n{} samples'.format(ref, len(mat_feature)))
plt.show()

### construct graph ###
import igraph as ig
import leidenalg as la

c_thresh = 0.25
i_vec, j_vec = np.triu_indices(len(mat_feature), k=1)
c_mask = c_vec>=c_thresh

sel_i_vec = i_vec[c_mask]
sel_j_vec = j_vec[c_mask]
sel_c_vec = c_vec[c_mask]
g = ig.Graph()
g.add_vertices(len(mat_feature))
g.add_edges(zip(sel_i_vec, sel_j_vec))
g.es['weight'] = sel_c_vec
partition = la.find_partition(g, la.ModularityVertexPartition, weights='weight')

ig.plot(
    g,
    layout=g.layout_auto(),
    target='/home/adrian/Data/DRACH/TAATCCCAGCACTTTGGGAGG/cluster.pdf'
)