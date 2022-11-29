import os
import pickle
from Bio.pairwise2 import align, format_alignment
import numpy as np
from collections import Counter
from tqdm import tqdm
from sklearn.metrics import pairwise_distances as sk_pdist
import matplotlib.pyplot as plt

HOME = os.environ['HOME']

motif = 'GGACT'
vocab = { 1:"A", 2:"C", 3:"G", 4:"T" }

def collect_motif_feature(loss_pred_label_feature):
    motifs = []
    features = []
    losses = []
    for event in tqdm(loss_pred_label_feature.values()):
        loss = event['loss']
        pred_label = event['pred_label']
        gt_label = event['gt_label']
        feature = event['feature']

        alignments = align.globalxx(gt_label, pred_label)
        aligned = alignments[0]

        motif_start_pos = aligned.seqA.find(motif[::-1])
        if motif_start_pos < 0:
            continue
        ref_mod_pos = motif_start_pos + 2
        # if (aligned.seqA[ref_mod_pos] != 'A') or (pred_label[feat_mod_pos] != 'A'):
        if (aligned.seqA[ref_mod_pos] != 'A') or (aligned.seqB[ref_mod_pos] != 'A'):
            continue

        feat_mod_pos = ref_mod_pos - Counter(aligned.seqB[:ref_mod_pos+1])['-']

        mod_motif = pred_label[(feat_mod_pos - 2):(feat_mod_pos + 3)]
        if len(mod_motif)<5:
            continue
        mod_feature = feature[feat_mod_pos]
        features.append(mod_feature)
        motifs.append(mod_motif)
        losses.append(loss)
    features = np.vstack(features)
    return motifs, features, losses

img_out = os.path.join(HOME, 'img_out/DRACH')
os.makedirs(img_out, exist_ok=True)

feature_file_ivt = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293_IVT_2/DRACH/{}_loss_pred_gt_feature.pkl'.format(motif))
feature_file_wt = os.path.join(HOME, 'Data/Isabel_IVT_Nanopore/HEK293A_wildtype/DRACH/{}_loss_pred_gt_feature.pkl'.format(motif))

with open(feature_file_ivt, 'rb') as f_ivt:
    collection_ivt = pickle.load(f_ivt)
motifs_ivt, features_ivt, losses_ivt = collect_motif_feature(collection_ivt)

with open(feature_file_wt, 'rb') as f_wt:
    collection_wt = pickle.load(f_wt)
motifs_wt, features_wt, losses_wt = collect_motif_feature(collection_wt)

### histograms of losses ###
# min_edge = 0
# max_edge = 7
# bin_width = 0.1
# num_bins = int((max_edge - min_edge) / bin_width)
# bin_edges = np.linspace(min_edge, max_edge, num_bins+1)
# bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
#
# hist_loss_wt, _ = np.histogram(losses_wt, bins=bin_edges)
# pdf_loss_wt = hist_loss_wt / np.sum(hist_loss_wt)
#
# hist_loss_ivt, _ = np.histogram(losses_ivt, bins=bin_edges)
# pdf_loss_ivt = hist_loss_ivt / np.sum(hist_loss_ivt)
#
# plt.figure(figsize=(10, 8))
# plt.plot(bin_centers, pdf_loss_ivt, label='HEK293 IVT 2')
# plt.plot(bin_centers, pdf_loss_wt, label='HEK293A WT')
# plt.xlabel('CTC loss')
# plt.ylabel('Density')
# plt.legend(loc='upper right')
# plt.savefig(os.path.join(img_out, 'hist_ctc_loss_wt_vs_ivt.png'), bbox_inches='tight')
# plt.close()

### calculate pairwise distance ###
# too slow!!!
# from scipy.spatial.distance import pdist
# c_vec = 1.0 - pdist(mat_feature, 'cosine')

ds_label = np.array(['IVT' for i in range(features_ivt.shape[0])] + ['WT' for i in range(features_wt.shape[0])])
mat_feature = np.vstack((features_ivt, features_wt))
c_mat = 1.0 - sk_pdist(X=mat_feature, metric='cosine', n_jobs=-1)
i_vec, j_vec = np.triu_indices(len(mat_feature), k=1)
c_vec = c_mat[i_vec, j_vec]

### plot histogram ###
c_thresh = 0.45

plt.figure(figsize=(10, 6))
plt.hist(c_vec, bins=100, range=[-1, 1])
plt.xlabel('Cosine similarity')
plt.ylabel('Counts')
plt.axvline(c_thresh, color='r')
plt.title('{}\n{} samples'.format(motif, len(mat_feature)))
plt.savefig(os.path.join(img_out, '{}_hist_similarity.png'.format(motif)), bbox_inches='tight')
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

c_thresh = 0.45
i_vec, j_vec = np.triu_indices(len(mat_feature), k=1)
c_mask = c_vec>=c_thresh

sel_i_vec = i_vec[c_mask]
sel_j_vec = j_vec[c_mask]
sel_c_vec = c_vec[c_mask]
g = ig.Graph()
g.add_vertices(len(mat_feature))
g.vs['ds_label'] = ds_label
g.add_edges(zip(sel_i_vec, sel_j_vec))
g.es['weight'] = sel_c_vec
partition = la.find_partition(g, la.ModularityVertexPartition, weights='weight')
membership = np.array(partition.membership)
modularity = partition.modularity

# num_clusters = len(np.unique(membership))
# palette = ig.RainbowPalette(n=num_clusters)
# for cluster_ind in range(num_clusters):
#     event_ind = np.where(membership==cluster_ind)[0]
#     g.vs[event_ind]["color"] = cluster_ind
#     cluster_edges = g.es.select(_within=g.vs[event_ind])
#     cluster_edges["color"] = cluster_ind
#
# fig = plt.figure(figsize=(10, 5))
# ax1 = fig.add_subplot(1, 2, 1)
# ig.plot(partition, target=ax1)

ig.plot(partition, os.path.join(img_out, '{}_clusters.png'.format(motif)))

# ig.plot(
#     g,
#     layout=g.layout_fruchterman_reingold(),
#     bbox=(100, 100),
#     palette=palette,
#     edge_width=1,
#     target=ax1,
#     vertex_size=0.1,
# )
plt.title('Modularity {:.2f}'.format(modularity))

for cluster_ind in range(num_clusters):
    event_ind = np.where(membership==cluster_ind)[0]
    g.vs[event_ind]["color"] = cluster_ind
    cluster_edges = g.es.select(_within=g.vs[event_ind])
    cluster_edges["color"] = cluster_ind

ax2 = fig.add_subplot(1, 2, 2)
ig.plot(
    g,
    layout=g.layout_fruchterman_reingold(),
    bbox=(100, 100),
    edge_width=1,
    target=ax1,
    vertex_size=0.1,
)
plt.title('Modularity {:.2f}'.format(modularity))


plt.savefig(os.path.join(img_out, '{}_clusters.png'.format(motif)), bbox_inches='tight')

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