"""
I noticed I was using the same few lines of code over and over to generate plots,
so I thought I would write a few functions to generate and save these plots.
"""

import matplotlib.pyplot as plt


#positions is a list of position tuples: [(x1,y1), (x2,y2), ... ]
def simple_graph_1d(positions, name = None, lw = 0.5):
    plt.plot(positions, linewidth = lw)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    if name:
        plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
    else: plt.show()

def simple_graph_2d(positions, name = None, lw = 0.5):
    x, y = zip(*positions)
    plt.plot(x, y, linewidth = lw)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    if name:
        plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
    else: plt.show()