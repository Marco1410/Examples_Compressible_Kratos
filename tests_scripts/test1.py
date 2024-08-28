from matplotlib import pyplot as plt
from scipy.stats import norm
import numpy as np

def calculate_noisy_line(grad, intcpt, mean_1, scale_1, mean_2, scale_2):
    # xs in range -150 to 150
    xs = np.linspace(-150, 150, 5001)
    # non-noisy line
    ys = grad * xs + intcpt
    # add noise from two Gaussians
    noise = norm.pdf(xs, mean_1, scale_1) + norm.pdf(xs, mean_2, scale_2)
    return xs, ys + noise

# Calculate and plot a noisy line
def plot_noisy_line(axis, grad, intcpt, mean_1, scale_1, mean_2, scale_2):
    xs, ys = calculate_noisy_line(grad, intcpt, mean_1, scale_1, mean_2, scale_2)
    axis.plot(xs, ys)

# Plot the two lines we're using in this notebook (fixed noise)
def plot_noisy_lines(axis):
    plot_noisy_line(axis, 0.02, 1, -4, 1, 3, 2)
    plot_noisy_line(axis, 0.019, 0.99, -3.9, 1.1, 2.7, 1.75)

# Just plot the data as is
def plot_original(axis):
    plot_noisy_lines(axis)
    axis.set_xlim(-150, 150)
    axis.set_title('Original')

def plot_delta(axis):
    xs, y1s = calculate_noisy_line(0.02, 1, -4, 1, 3, 2)
    _, y2s = calculate_noisy_line(0.019, 0.99, -3.9, 1.1, 2.7, 1.75)
    deltas = np.abs(y2s - y1s)
    axis.plot(xs, deltas,
              color=plt.rcParams['axes.prop_cycle'].by_key()['color'][2])

def plot_different_data(axis):
    plot_original(axis)
    axis.set_title('Different Data')

    # Inset axes showing the absolute difference between the two lines
    axin1 = axis.inset_axes([0.2, 0.6, 0.2, 0.2])
    plot_delta(axin1)
    axin1.set_xlim(0, 0.01)
    axin1.set_xticks([-15, 0, 15])
    axin1.set_xticklabels(['-1 ', '0', '1'])
    axin1.set_ylim(0, 0.01)
    axin1.set_yticks([0, 0.15])
    axin1.set_yticklabels(['0 ', '0.15'])
    axin1.set_title('Abs. Diff.', size=13, pad=8)

    # Create the usual zoomed in inset axes in the top left corner
    # axin2 = axis.inset_axes([0.04, 0.6, 0.36, 0.36])
    # plot_noisy_lines(axin2)
    # axin2.set_xlim(-10, 12)
    # axin2.set_ylim(0.75, 1.4)
    # axin2.set_xticks([])
    # axin2.set_yticks([])
    axis.indicate_inset_zoom(axin1)

fig, axis = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
plot_different_data(axis)
plt.show()