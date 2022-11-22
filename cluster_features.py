import os
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
img_out = '/home/adrian/img_out/DRACH'
os.makedirs(img_out, exist_ok=True)

vocab = { 1:"A", 2:"C", 3:"G", 4:"T" }

with open(file, 'rb') as f:
    loss_pred_label_feature = pickle.load(f)

motifs = []
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
        # print(aligned.seqA[:ref_mod_pos+1])
        # print(pred_label[:pred_mod_pos+1])
        # print('\n')

        mod_feature = feature[pred_mod_pos]
        mat_feature.append(mod_feature)

        motifs.append(pred_label[(pred_mod_pos-2):(pred_mod_pos+3)])

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

c_thresh = 0.20

plt.figure(figsize=(10, 6))
plt.hist(c_vec, bins=50, range=[-1, 1])
plt.xlabel('Cosine similarity')
plt.ylabel('Counts')
plt.axvline(c_thresh, color='r')
plt.title('{}\n{} samples'.format(ref, len(mat_feature)))
plt.savefig(os.path.join(img_out, 'hist_similarity.png'), bbox_inches='tight')
# plt.show()

### connected components ###
# from scipy.sparse import coo_matrix
# from scipy.sparse.csgraph import connected_components
#
# i_vec, j_vec = np.triu_indices(len(mat_feature), k=1)
# # c_mat = coo_matrix((c_vec, (i_vec, j_vec)), shape=(len(mat_feature), len(mat_feature)))
# #
# # plt.figure(figsize=(8, 8))
# # plt.imshow(c_mat.toarray())
# # plt.show()
#
# thresh_ncpts_G = {}
# for c_thresh in np.arange(0, 1, 0.01):
#     c_mask = c_vec >= c_thresh
#     sel_i_vec = i_vec[c_mask]
#     sel_j_vec = j_vec[c_mask]
#     sel_c_vec = c_vec[c_mask]
#     c_mat = coo_matrix((sel_c_vec, (sel_i_vec, sel_j_vec)), shape=(len(mat_feature), len(mat_feature)))
#     n_cpts, membership = connected_components(c_mat, directed=False)
#     G = Counter(membership).most_common(1)[0][1]
#     thresh_ncpts_G[c_thresh] = (n_cpts, G)
#
# plt.figure(figsize=(10, 6))
# plt.plot(thresh_ncpts_G.keys(), [val[0] for val in thresh_ncpts_G.values()], label='N')
# plt.plot(thresh_ncpts_G.keys(), [val[1] for val in thresh_ncpts_G.values()], label='G')
# plt.xlabel('C threshold')
# # plt.ylabel('Connected components')
# plt.legend(loc='upper left')
# plt.savefig('/home/adrian/Data/DRACH/TAATCCCAGCACTTTGGGAGG/thresh_conn_cpts.png', bbox_inches='tight')
# plt.show()

### fix c thresh ###
# c_thresh = 0.80
# c_mask = c_vec >= c_thresh
# sel_i_vec = i_vec[c_mask]
# sel_j_vec = j_vec[c_mask]
# sel_c_vec = c_vec[c_mask]
# c_mat = coo_matrix((sel_c_vec, (sel_i_vec, sel_j_vec)), shape=(len(mat_feature), len(mat_feature)))
# n_cpts, membership = connected_components(c_mat, directed=False)


### construct graph ###
import igraph as ig
import leidenalg as la

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
membership = np.array(partition.membership)
modularity = partition.modularity

num_clusters = len(np.unique(membership))
palette = ig.RainbowPalette(n=num_clusters)
for cluster_ind in range(num_clusters):
    event_ind = np.where(membership==cluster_ind)[0]
    g.vs[event_ind]["color"] = cluster_ind
    cluster_edges = g.es.select(_within=g.vs[event_ind])
    cluster_edges["color"] = cluster_ind

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ig.plot(
    g,
    layout = g.layout_fruchterman_reingold(),
    bbox = (100, 100),
    palette=palette,
    edge_width=1,
    target=ax,
    vertex_size=0.1,
)
plt.title('Modularity {:.2f}'.format(modularity))
plt.savefig(os.path.join(img_out, 'clusters.png'), bbox_inches='tight')

### outlier ###
counts = Counter(membership).most_common()
# outlier_indices = np.where(membership==counts[-1][0])[0]
# outlier_motifs = [motifs[i] for i in outlier_indices]
# inlier_indices = [i for i in range(len(motifs)) if i not in outlier_indices]
# inlier_motifs = [motifs[i] for i in inlier_indices]

print('Thresh {}:'.format(c_thresh))
for cluster_ind, cluster_size in counts:
    print('Cluster {}, {} reads'.format(cluster_ind, cluster_size))
    indices = np.where(membership==cluster_ind)[0]
    index_motifs = np.array(motifs)[indices]
    for k, v in Counter(index_motifs).most_common():
        print('{}: {}'.format(k[::-1], v))