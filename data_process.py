import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
import argparse
import pandas as pd
from scipy.spatial import ConvexHull
import pickle as pkl

def segment_no_events(data, no_event_val=-1):
    no_events_regions = [r.astype('int') for r in measure.find_contours(data, no_event_val)]
    no_events_regions_sizes = np.array([r.shape[0] for r in no_events_regions])
    sorted_regions_idx = np.argsort(no_events_regions_sizes)
    inner_region = no_events_regions[sorted_regions_idx[-2]]
    outer_region = no_events_regions[sorted_regions_idx[-1]] # the outer region is the largest
    center = inner_region.mean(axis=0)
    return inner_region, outer_region, center

parser = argparse.ArgumentParser(
                    prog='Thin film Data Analyzer')

parser.add_argument('--field', type=float, default=0.13) 
parser.add_argument('--basename', default="raw_data/file_field")
parser.add_argument('--n_th', type=int, default=2**9) 

args = parser.parse_args()


filename = args.basename + ("%1.2f" % args.field) + ".csv"
data = np.loadtxt(filename, delimiter=",").astype('int')

n_rows, n_cols = data.shape[0], data.shape[1]

inner_region, outer_region, center = segment_no_events(data,  no_event_val=-1)
ii, jj = np.meshgrid(np.arange(0,n_rows), np.arange(0,n_cols), indexing="ij")



data_df = pd.DataFrame()
data_df['row'] = ii.flatten()
data_df['col'] = jj.flatten()
data_df['frame'] = data.flatten()
data_df['r'] = np.linalg.norm(data_df[['row', 'col']] - center[None,:] )
data_df['th'] = np.arctan2(data_df['row'] - center[0],data_df['col'] - center[1]  )
data_df = data_df[data_df.frame != -1]
#data_df.drop(data_df.frame == -1, axis=0, inplace=True)
data_df.reset_index(drop=True, inplace=True)

xs = {}
ys = {}
for frame, rows in data_df.groupby('frame'):
    xs[frame] = rows.row.values
    ys[frame] = rows.col.values
with open(args.basename +  ("%1.2f" % args.field) + "_events.pkl",'wb') as file:
    pkl.dump([xs, ys], file)
    
    
    



exit()

# interface reconstruction


pts = np.vstack([ii.flatten(), jj.flatten()])


outer_hull = ConvexHull(outer_region)
outer_mask = np.ones((n_rows, n_cols), dtype=bool)
for eq in outer_hull.equations:
    eq_mask = np.dot(eq[:-1], pts) + eq[-1] <= 1e-12
    outer_mask[ii.flatten(),jj.flatten()] *= eq_mask

inner_hull = ConvexHull(inner_region)
inner_mask = np.ones((n_rows, n_cols), dtype=bool)
for eq in inner_hull.equations:
    eq_mask = np.dot(eq[:-1], pts) + eq[-1] <= 1e-12
    inner_mask[ii.flatten(),jj.flatten()] *= eq_mask

no_crown_mask = ~(outer_mask*(1-inner_mask))
















n_thetas = args.n_th
    