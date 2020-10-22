import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.font_manager
from skimage import color

import warnings


if 'Fira' not in matplotlib.font_manager.findfont('Fira Sans'):
    warnings.warn('Fira Sans font not installed, or matplotlib font cache out '
                  'of date.')

col_width = 3.405  # in, size of figure width=\columnwidth
true_col_width = 3.463  # in, size of example column of text
text_width = 6.998  # in, size of figure* width=\textwidth
true_text_width = 7.024  # in, size of two example columns together
content_width = 7.518  # in, including line numbers, etc
page_width = 8.5  # 8.5 x 11in format
golden_ratio = (1 + np.sqrt(5))/2


def use_pnas_style():
    # PNAS declares 9pt default text but uses ~8pt, with ~6pt captions
    mpl.rcParams['font.size'] = 8.0
    # however, Sean kept requiresting larger text, so now we're up to
    # 9.6pt font (large), and ~11.52pt titles (x-large)
    mpl.rcParams['xtick.labelsize'] = 'large'  # default : medium
    mpl.rcParams['ytick.labelsize'] = 'large'  #
    mpl.rcParams['axes.titlesize'] = 'x-large'  # default : large,
    mpl.rcParams['axes.labelsize'] = 'large'  # default : medium
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['font.serif'] = 'Fira Sans'
    mpl.rcParams['font.family'] = 'sans-serif'

# Generate color maps, ScalarMappables, tickers, locators, etc.

vmax_long = 5.7
vmax_light = 5.0
vmax_sim = 7.5
total_colors = 5000

two_colors = {'blue': [63.5/100.0, 82.4/100.0, 100.0/100.0],
              'purple': [80.4/100.0, 70.6/100.0, 85.9/100.0]}
grey = np.array([100/255, 100/255, 100/255])
three_colors = [[ 56/255   , 125/255   , 156/255  ],
                [ 63/255   , 138/255   ,  92/255  ],
                [163/255   , 85/255    , 153/255  ]]
# #                 [166/255   , 58/255    , 133/255  ]]
# need to add a dim, since skimage expects >=2d images
lab_colors = color.rgb2lab(np.array([three_colors]), illuminant='D65', observer='2')[0]
lab_sim_colors = color.rgb2lab(np.array([list(two_colors.values())]),
                               illuminant='D65',
                               observer='2')[0]

def bruno_div_map(n, left_lab, right_lab, midpoint_l, target_l=None):
    l, a, b = 0, 1, 2
    left_l = left_lab[l] if target_l is None else target_l
    right_l = right_lab[l] if target_l is None else target_l
    t = np.linspace(0, 1, n)
    lab_a = np.interp(t, [0, 1], [left_lab[a], right_lab[a]])
    lab_b = np.interp(t, [0, 1], [left_lab[b], right_lab[b]])
    lab_l = np.interp(t, [0, 1/2, 1], [left_l, midpoint_l, right_l])
    return list(zip(lab_l, lab_a, lab_b))

# artificially place the sawtooth peak at T=3.5 to highlight "mid-prophase",
# but also cause the "0.5" extra makes teh post-discretized colors look nice
num_grey = int(((0 - (-1))/(vmax_long - (-1))) * total_colors) - 1
num_left = int(((3.5 - 0)/(vmax_long - (-1))) * total_colors)
num_right = int(((vmax_long - 3.5)/(vmax_long - (-1))) * total_colors)
num_sim = 15  # whatever, man, matplotlib smoothes

blended_divs = np.concatenate((
    bruno_div_map(num_left, lab_colors[0], lab_colors[1], midpoint_l=80, target_l=50),
    bruno_div_map(num_right, lab_colors[1], lab_colors[2], midpoint_l=80, target_l=50)[1:]
))
# need to add a dim, since skimage expects >=2d images
bgr_long = color.lab2rgb(np.array([blended_divs]), illuminant='D65', observer='2')[0]
bgr_long = np.concatenate((np.tile(grey, (num_grey, 1)), bgr_long))

num_grey = int(((0 - (-1))/(vmax_light - (-1))) * total_colors) - 1
num_left = int(((3 - 0)/(vmax_light - (-1))) * total_colors)
num_right = int(((vmax_light - 3)/(vmax_light - (-1))) * total_colors)

blended_divs_light = np.concatenate((
    bruno_div_map(num_left, lab_colors[0], lab_colors[1], midpoint_l=80, target_l=50),
    bruno_div_map(num_right, lab_colors[1], lab_colors[2], midpoint_l=80, target_l=50)[1:]
))
# need to add a dim, since skimage expects >=2d images
bgr_light = color.lab2rgb(np.array([blended_divs_light]), illuminant='D65', observer='2')[0]
bgr_light = np.concatenate((np.tile(grey, (num_grey, 1)), bgr_light))

# Co-authors want a separate, linear cmap
blended_sim_divs = bruno_div_map(num_sim, lab_sim_colors[0], lab_sim_colors[1], midpoint_l=80, target_l=70)
bgr_sim = color.lab2rgb(np.array([blended_sim_divs]), illuminant='D65', observer='2')[0]
# sim_cmap = mpl.cm.winter_r
sim_cmap = mpl.colors.ListedColormap(bgr_sim)


cmap_long = mpl.colors.ListedColormap(bgr_long)
cmap_light = mpl.colors.ListedColormap(bgr_light)
# sim_cmap = mpl.colors.ListedColormap(bgr_dark)
# sim_cmap = mpl.colors.ListedColormap(bgr_sim)


# discretize the continuous norm to integers in the relevant range
long_labels = np.arange(-1, np.ceil(vmax_long) + 1)
long_cbounds = np.arange(-1.5, np.ceil(vmax_long) + 1.5)
# continuous norm, to be replaced once we extract out the discrete colors
long_cnorm_continuous = mpl.colors.Normalize(vmin=-1, vmax=vmax_long)
long_colors = cmap_long(long_cnorm_continuous(long_labels))
long_cnorm = mpl.colors.BoundaryNorm(long_cbounds, len(long_labels))
long_listed_cmap = mpl.colors.ListedColormap(long_colors)
long_sm = mpl.cm.ScalarMappable(norm=long_cnorm, cmap=long_listed_cmap)
long_sm.set_array([])
long_locator = ticker.MultipleLocator(1)
long_formatter = ticker.FuncFormatter(lambda i, _: f'$T_{{{int(i)}}}$' if i >= 0 else '$G_0$')
long_sm_continuous = mpl.cm.ScalarMappable(norm=long_cnorm_continuous, cmap=cmap_long)
long_sm_continuous.set_array([])

sim_labels = np.arange(1, vmax_sim + 1)
sim_cbounds = np.arange(0.5, vmax_sim + 1.5)
# sim_cbounds[0] = 1
sim_cnorm_continuous = mpl.colors.Normalize(vmin=1, vmax=vmax_sim)
sim_colors = sim_cmap(sim_cnorm_continuous(sim_labels))
sim_cnorm = mpl.colors.BoundaryNorm(sim_cbounds, len(sim_labels))
sim_listed_cmap = mpl.colors.ListedColormap(sim_colors)
sim_locator = ticker.MultipleLocator(1)
sim_sm = mpl.cm.ScalarMappable(norm=sim_cnorm, cmap=sim_listed_cmap)
sim_sm.set_array([])
sim_sm_continuous = mpl.cm.ScalarMappable(norm=sim_cnorm_continuous, cmap=sim_cmap)
sim_sm_continuous.set_array([])

# Ignore for now, see what breaks.
# cnorm_light = mpl.colors.Normalize(vmin=-1, vmax=vmax_light)
# # scalar mappables are just norm/cmap pairs
# meiosis_sm_light = mpl.cm.ScalarMappable(norm=cnorm_light, cmap=cmap_light)
# meiosis_sm_light.set_array([])

cycling_cmap = plt.get_cmap("tab10")
