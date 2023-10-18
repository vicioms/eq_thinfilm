import numpy as np
import matplotlib.pyplot as plt
import argparse
import os



parser = argparse.ArgumentParser(
                    prog='Thin film Data Analyzer')

parser.add_argument('--field', type=float, default=0.13) 
parser.add_argument('--basename', default="data/file_field")

args = parser.parse_args()

intf_filename = args.basename +  ("%1.2f" % args.field) + "_intf.npz"


interface = np.load(intf_filename)[0]
print(interface)