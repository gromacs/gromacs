import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def plot(res, legend=None, title=None,
         xlabel=None, ylabel=None, xlim=None, ylim=None,
         inv_x=False, inv_y=False,
         filename=None, screen=True):

    font = {'family': 'serif',
            'weight': 'normal',
            'size': 16}
    mpl.rc('font', **font)

    plt.ioff()
    fig, ax = plt.subplots()
    xmin = float('inf')
    xmax = float('-inf')
    for r in res:
        x = r['x']
        y = r['y']
        if xlim is not None:
            x = x[(r['x'] >= xlim[0]) & (r['x'] <= xlim[1])]
            y = y[(r['x'] >= xlim[0]) & (r['x'] <= xlim[1])]
        if 'y_err' in r:
            dy = r['y_err']
            if xlim is not None:
                dy = dy[(r['x'] >= xlim[0]) & (r['x'] <= xlim[1])]
            ax.errorbar(x, y, label=r['name'], yerr=dy)
        else:
            ax.plot(x, y, label=r['name'])
        if xlim is not None:
            xmin = min(np.min(r['x']), xmin)
            xmax = max(np.max(r['x']), xmax)
        else:
            xmin = min(np.min(r['x']), xmin)
            xmax = max(np.max(r['x']), xmax)

    if legend is not None:
        ax.legend(loc=legend)
    box = ax.get_position()
    if title is not None:
        ax.set_title(title, y=1.05)
        box = box.from_bounds(box.x0, box.y0, box.width, box.height * 0.95)
    if xlabel is not None:
        ax.set_xlabel(xlabel, labelpad=5)
        box = box.from_bounds(box.x0, box.y0 + 0.05 * box.height, box.width, box.height * 0.95)
    if ylabel is not None:
        ax.set_ylabel(ylabel, labelpad=10)
        box = box.from_bounds(box.x0 + 0.05 * box.width, box.y0, box.width * 0.95, box.height)
    ax.set_position([box.x0, box.y0, box.width, box.height])
    ax.axis('auto')
    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim([xmin, xmax])
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    if inv_x:
        ax.invert_xaxis()
    if inv_y:
        ax.invert_yaxis()

    ax.ticklabel_format(style='sci', axis='x', scilimits=(-3, 4))
    ax.xaxis.major.formatter._useMathText = True

    if filename is not None:
        fig.savefig(filename + '.pdf', dpi=300)
    if screen:
        fig.show()
