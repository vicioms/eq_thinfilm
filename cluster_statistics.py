import numpy as np
import networkx as nx
import pickle as pkl
import argparse
import scipy.sparse as sp
import matplotlib.pyplot as plt

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
sparse_net[sparse_net > args.ds] = 0
net = nx.Graph(sparse_net)
clusters = []
comps = nx.connected_components(net)
for comp in comps:
    clusters.append(np.sort(list(comp)))
print("Number of clusters:", len(clusters))
cl_sizes = []
for cluster in clusters:
    if(len(cluster) == 1 and cluster[0] not in xs):
        continue
    x_c = []
    y_c = []
    for c in cluster:
        if(c in xs):
            x_c.append(xs[c])
            y_c.append(ys[c])
    x_c = np.concatenate(x_c)
    y_c = np.concatenate(y_c)
    cl_sizes.append(len(x_c))
cl_sizes = np.array(cl_sizes)
plt.plot(cl_sizes)
plt.yscale('log')
plt.show()
hist, edges = np.histogram(cl_sizes, bins=np.logspace(0,1.1*np.log10(np.max(cl_sizes)),30),density=True )
plt.scatter(edges[:-1],hist)
plt.plot(edges, 0.1*edges**(-1.26))
plt.plot(edges, 0.1*edges**(-1.11))
plt.xscale('log')
plt.yscale('log')
plt.show()