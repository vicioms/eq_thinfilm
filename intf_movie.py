import numpy as np
import matplotlib.pyplot as plt
import argparse
import os



parser = argparse.ArgumentParser(
                    prog='Thin film Data Analyzer')

parser.add_argument('--field', type=float, default=0.13) 
parser.add_argument('--basename', default="data/file_field")

args = parser.parse_args()

intf_filename = args.basename +  ("%1.2f" % args.field) + "_intf.npy"


interface = np.load(intf_filename)


import matplotlib.pyplot as plt
import numpy as np

import matplotlib.animation as animation

fig, ax = plt.subplots()
ax.set_ylim(-25,25)
x = np.linspace(-np.pi, np.pi, interface.shape[1])
line, = ax.plot(x, interface[0]-interface[0].mean(), color='black')


def animate(i):
    line.set_ydata(interface[i,:]-interface[i,:].mean())  # update the data.
    return line,


ani = animation.FuncAnimation(
    fig, animate, interval=2, blit=True, save_count=50)

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# writer = animation.FFMpegWriter(
#     fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)

plt.show()
