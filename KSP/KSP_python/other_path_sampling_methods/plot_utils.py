import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

def plot_hexbin_with_max(
    x, y, C,
    gridsize=40,
    mincnt=1,
    vmin=0.0, vmax=0.3,
    cmap_name='inferno_r',
    max_marker_radius=200,
    annotate_max=True,
    axis_off=True,
    fig=None, ax=None,
    cbar_label = None
):
    """
    Hexbin plot with values C, highlighting the cell with the maximum C.

    Parameters
    ----------
    x, y : array-like
        Coordinates (same length).
    C : array-like
        Values to aggregate per hex (same length as x,y).
    gridsize : int or (int, int)
        Hexbin resolution.
    mincnt : int
        Only display hexes with at least this many points.
    vmin, vmax : float
        Color normalization limits.
    cmap_name : str
        Matplotlib colormap name (e.g., 'inferno_r', 'RdBu_r').
    max_marker_radius : float
        Radius of the circle highlighting the max-valued hex (in data coords).
    annotate_max : bool
        If True, add legend entry showing the max value.
    axis_off : bool
        If True, turn off axes frame/ticks.
    fig, ax : matplotlib Figure/Axes
        Optional existing fig/ax to draw on.
    cbar_label : str
        Label for the colorbar.

    Returns
    -------
    fig, ax, hb : (Figure, Axes, PolyCollection)
        The created (or passed-in) figure, axes, and hexbin artist.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    C = np.asarray(C)

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    if type(cmap_name) == str:
        cmap = plt.colormaps[cmap_name]
    else:
        cmap = cmap_name

    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    hb = ax.hexbin(
        x, y, C=C, gridsize=gridsize, cmap=cmap, norm=norm,
        mincnt=mincnt, alpha=0.9
    )

    # Find max-valued hex and outline it
    offs = hb.get_offsets()       # (Nhex, 2) centers of hexes
    vals = hb.get_array()         # (Nhex,) aggregated C values per hex
    print('max value:', np.max(vals), 'at position:', offs[np.argmax(vals)])
    if vals.size > 0 and np.isfinite(vals).any():
        max_ind = np.argmax(vals)
        max_val = np.max(vals)
        circle=plt.Circle(offs[max_ind], max_marker_radius, color = 'green',
                           fill=False, linewidth=3, 
                           label = 'max val = '+str(round(max_val, 3))) 

    # Colorbar
    cbar = plt.colorbar(hb, ax=ax, orientation='vertical')
    cbar.set_label(cbar_label, fontsize=30)
    cbar.ax.tick_params(labelsize=25)

    # Cosmetics
    if axis_off:
        ax.axis('off')
    ax.set_aspect('equal', adjustable='box')

    if vals.size > 0 and np.isfinite(vals).any():
        ax.add_patch(circle)

    if annotate_max:
        ax.legend(fontsize=16)

    plt.tight_layout()
    return fig, ax, hb


# -----------------------
# Example usage:
# fig, ax, hb = plot_hexbin_with_max(
#     x_pos, y_pos, np.array(sp_bwss_ls[0]),
#     gridsize=40, mincnt=1,
#     vmin=0.0, vmax=0.3,
#     cmap_name='inferno_r',
#     max_marker_radius=0.05
# )
# plt.show()
