import numpy as np
import pickle as pkl
import argparse
import scipy.sparse as sp

parser = argparse.ArgumentParser(
                    prog='Thin film Data Analyzer - Clusters')

parser.add_argument('--field', type=float, default=0.13) 
parser.add_argument('--dt', type=int, default=10)
parser.add_argument('--basename', default="raw_data/file_field")

args = parser.parse_args()
filename = args.basename + "%1.2f_events.pkl" % args.field
with open(filename, 'rb') as file:
    xs, ys = pkl.load(file)
frames = np.sort(list(xs.keys()))
max_frames = np.max(frames)+1

dt = args.dt

sparse_net = sp.dok_array(arg1=(max_frames, max_frames), dtype=int)
for counter, f1 in enumerate(frames):
    x1, y1 = xs[f1], ys[f1]
    for f2 in frames[counter+1:]:
        assert f2 > f1
        if(f2-f1 > dt):
            break
        x2, y2 = xs[f2], ys[f2]
        d_21 = np.min(np.abs(x1[None, :]-x2[:,None])+np.abs(y1[None, :]-y2[:,None]))
        sparse_net[f1, f2] = d_21
        sparse_net[f2, f1] = d_21
sp.save_npz(args.basename + "%1.2f_net_dt=%i" % (args.field, dt), sparse_net.tocsr())