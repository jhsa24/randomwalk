"""
I noticed I was using the same few lines of code over and over to generate plots,
so I thought I would write a few functions to generate and save these plots.
"""

import matplotlib.pyplot as plt


#positions is a list of position tuples: [(x1,y1), (x2,y2), ... ]
def simple_graph(positions, name = "simple_graph"):
    x, y = zip(*positions)
    plt.plot(x, y, linewidth = 0.3)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.savefig(name + ".png", dpi=300, bbox_inches='tight')