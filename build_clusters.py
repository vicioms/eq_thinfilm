import numpy as np
import networkx as nx
import pickle as pkl
import argparse

parser = argparse.ArgumentParser(
                    prog='Thin film Data Analyzer - Clusters')

parser.add_argument('--field', type=float, default=0.13) 
parser.add_argument('--ds', type=int, default=3)
parser.add_argument('--dt', type=int, default=10)
parser.add_argument('--basename', default="raw_data/file_field")

args = parser.parse_args()
filename = args.basename + "%1.2f_events.pkl" % args.field
with open(args.basename % args.field, 'rb') as file:
    xs, ys = pkl.load(file)
frames = list(xs.keys())


print(frames)
for f1 in frames:
    for f2 in frames:
        