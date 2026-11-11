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
import matplotlib.cm as cm

fig, ax = plt.subplots()
#ax.set_ylim(-5,5)
ax.set_ylim(50, 300)
x = np.linspace(-np.pi, np.pi, interface.shape[1])
line2, = ax.plot(x, interface[1], color='red')
line, = ax.plot(x, interface[0], color='black')
cmap = cm.get_cmap('jet')
max_frames = 50000

def animate(i):
    line2.set_ydata(interface[i+1,:])  # update the data.
    line.set_ydata(interface[max(i-100,0),:])
    #line.set_color(cmap(i/max_frames))
    return line,line2


ani = animation.FuncAnimation(
    fig, animate, interval=1, blit=True, frames=max_frames)

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
