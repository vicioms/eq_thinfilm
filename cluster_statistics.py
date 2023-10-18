import numpy as np
import networkx as nx
import pickle as pkl
import argparse
import scipy.sparse as sp
import matplotlib.pyplot as plt
from numba import njit
import pandas as pd


def cluster_length(x, y):
    c = np.zeros((2,2))
    c[0,0] = np.var(x)
    c[0,1] = np.mean(x*y)-np.mean(x)*np.mean(y)
    c[1,0] = c[0,1]
    c[1,1] = np.var(y)
    w = np.linalg.eigvalsh(c)
    return np.sqrt(w[1])

parser = argparse.ArgumentParser(
                    prog='Thin film Data Analyzer - Clusters')

parser.add_argument('--field', type=float, default=0.13) 
parser.add_argument('--dt', type=int, default=10)
parser.add_argument('--ds', type=int, default=2)
parser.add_argument('--basename', default="data/file_field")

args = parser.parse_args()
sparse_net_filename = args.basename + "%1.2f_net_dt=%i.npz" % (args.field, args.dt)
filename = args.basename + "%1.2f_events.pkl" % args.field
with open(filename, 'rb') as file:
    xs, ys = pkl.load(file)

sparse_net  = sp.load_npz(sparse_net_filename)
net = nx.Graph(sparse_net)
# filter out all edges
long_edges = list(filter(lambda e: e[2] > args.ds, (e for e in net.edges.data('weight'))))
le_ids = list(e[:2] for e in long_edges)
# remove filtered edges from graph
net.remove_edges_from(le_ids)
clusters = []
comps = nx.connected_components(net)
for comp in comps:
    clusters.append(np.sort(list(comp)))
print("Number of clusters:", len(clusters))
cl_sizes = []
cl_growth = []


for cluster in clusters:
    x_c = []
    y_c = []
    linear_stats = []
    total_pixels = 0
    for c in cluster:
        if(c in xs):
            total_pixels += len(xs[c])
            x_c.append(xs[c])
            y_c.append(ys[c])
            linear_stats.append((c, cluster_length(np.concatenate(x_c), np.concatenate(y_c)), total_pixels))
    if(len(x_c)==0):
        continue
    linear_stats = np.array(linear_stats)
    linear_stats[:, 0] -= linear_stats[0,0]
    cl_growth.append(linear_stats)
    x_c = np.concatenate(x_c)
    y_c = np.concatenate(y_c)
    cl_sizes.append(len(x_c))
cl_sizes = np.array(cl_sizes)
cl_growth = np.concatenate(cl_growth, axis=0)
df_growth = pd.DataFrame()
df_growth['timeIdx'] = cl_growth[:,0]
df_growth['linear'] = cl_growth[:,1]
df_growth['area'] = cl_growth[:,2]
df_growth['area_sq'] = cl_growth[:,2]**2
groups_mean = df_growth.groupby('timeIdx').mean()
plt.plot(np.arange(1, len(groups_mean)+1), groups_mean['area_sq']/ groups_mean['area'])
plt.xscale('log')
plt.yscale('log')
plt.show()

#hist, edges = np.histogram(cl_sizes, density=True, bins=np.logspace(0.2,1.1*np.log10(np.max(cl_sizes)),25) )
#plt.scatter(edges[:-1],hist)
#plt.plot(edges, 0.5*edges**(-1.26), color='black', ls='dashed',lw=3)
#plt.plot(edges, 0.5*edges**(-1.11), color='gray', ls='dashed', lw=3)
#plt.xscale('log')
#plt.yscale('log')
#plt.show()