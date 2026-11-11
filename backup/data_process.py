import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
import argparse
import pandas as pd
from scipy.spatial import ConvexHull
import pickle as pkl
import os


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
parser.add_argument('--basename', default="data/file_field")
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



pts = np.vstack([ii.flatten(), jj.flatten()])


outer_hull = ConvexHull(outer_region)
outer_mask = np.ones((n_rows, n_cols), dtype=bool)
for eq in outer_hull.equations:
    eq_mask = np.dot(eq[:-1], pts) + eq[-1] <= 1e-12
    outer_mask[ii.flatten(),jj.flatten()] *= eq_mask
outer_mask  = ~outer_mask
inner_hull = ConvexHull(inner_region)
inner_mask = np.ones((n_rows, n_cols), dtype=bool)
for eq in inner_hull.equations:
    eq_mask = np.dot(eq[:-1], pts) + eq[-1] <= 1e-12
    inner_mask[ii.flatten(),jj.flatten()] *= eq_mask
    
    
intdict_filename = args.basename +  ("%1.2f" % args.field) + "_intdict.pkl"

if(os.path.exists(intdict_filename)):
    print("Interface dictionary already exists.")
    with open(intdict_filename,'rb') as file:
        interface_dict = pkl.load(file)
else:
    print("Interface dictionary does not exist. Starting... (it might take a while)")
    temp_inner = np.zeros_like(data)
    temp_inner[inner_mask] = 1
    event_mask = data > - 1
    interface_dict = {}
    for t in range(0, data.max()+1):
        if(t%500 == 0):
            print("frame: ", t)
        temp = np.copy(temp_inner)
        temp[event_mask*(data<=t)] = 1
        regions = measure.find_contours(temp)
        r_sizes = [len(r) for r in regions]
        largest_r = regions[np.argmax(r_sizes)].astype('int')
        interface_dict[t] = largest_r

    with open(intdict_filename,'wb') as file:
        pkl.dump(interface_dict, file)
        
        
theta_new = np.linspace(-np.pi, np.pi, args.n_th)
interface_data_interp = []
for t, h in interface_dict.items():
    theta = np.arctan2(h[:, 0] - center[0], h[:, 1]-center[1])
    rho = np.sqrt((h[:, 0] - center[0])**2+(h[:, 1]-center[1])**2)
    rho_new =  np.interp(theta_new, theta, rho, period=np.pi*2)
    interface_data_interp.append(rho_new)
interface_data_interp = np.array(interface_data_interp)

interface = [np.copy(interface_data_interp[0])]
for t in range(1, interface_data_interp.shape[0]):
    new_interface = np.copy(interface_data_interp[t])
    mask_receeded = new_interface-interface[-1] < 0
    new_interface[mask_receeded] = interface[-1][mask_receeded]
    interface.append(new_interface)
    #print("Fixed sites: ", mask_receeded.sum())
interface = np.array(interface)

intf_filename = args.basename +  ("%1.2f" % args.field) + "_intf.npy"
np.save(intf_filename, interface)
exit()


# interface reconstruction

pts = np.vstack([ii.flatten(), jj.flatten()])


outer_hull = ConvexHull(outer_region)
outer_mask = np.ones((n_rows, n_cols), dtype=bool)
for eq in outer_hull.equations:
    eq_mask = np.dot(eq[:-1], pts) + eq[-1] <= 1e-12
    outer_mask[ii.flatten(),jj.flatten()] *= eq_mask
outer_mask  = ~outer_mask
inner_hull = ConvexHull(inner_region)
inner_mask = np.ones((n_rows, n_cols), dtype=bool)
for eq in inner_hull.equations:
    eq_mask = np.dot(eq[:-1], pts) + eq[-1] <= 1e-12
    inner_mask[ii.flatten(),jj.flatten()] *= eq_mask


masked_data = np.copy(data).astype('float')
masked_data[inner_mask] = 0
masked_data[outer_mask] = np.infty

frames = np.sort(list(xs.keys()))
frames = frames[frames>=0]
non_bg_mask = data != -1 + inner_mask
for frame in frames:
    data_next = np.zeros_like(masked_data)
    mask = (masked_data <= frame)*non_bg_mask
    data_next[mask] = 1
    if(frame % 1000 == 0):
        plt.imshow(data_next)
        plt.show()
























n_thetas = args.n_th
    