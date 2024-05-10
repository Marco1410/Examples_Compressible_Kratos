import math
import numpy
import matplotlib
import matplotlib.pyplot as plt

def draw_loglog_slope(fig, ax, origin, width_inches, slope, inverted=False, color=None, polygon_kwargs=None, label=True, labelcolor=None, label_kwargs=None, zorder=None):
    """
    This function draws slopes or "convergence triangles" into loglog plots.
    @param fig: The figure
    @param ax: The axes object to draw to
    @param origin: The 2D origin (usually lower-left corner) coordinate of the triangle
    @param width_inches: The width in inches of the triangle
    @param slope: The slope of the triangle, i.e. order of convergence
    @param inverted: Whether to mirror the triangle around the origin, i.e. whether
        it indicates the slope towards the lower left instead of upper right (defaults to false)
    @param color: The color of the of the triangle edges (defaults to default color)
    @param polygon_kwargs: Additional kwargs to the Polygon draw call that creates the slope
    @param label: Whether to enable labeling the slope (defaults to true)
    @param labelcolor: The color of the slope labels (defaults to the edge color)
    @param label_kwargs: Additional kwargs to the Annotation draw call that creates the labels
    @param zorder: The z-order value of the triangle and labels, defaults to a high value
    """

    if polygon_kwargs is None:
        polygon_kwargs = {}
    if label_kwargs is None:
        label_kwargs = {}

    if color is not None:
        polygon_kwargs["color"] = color
    if "linewidth" not in polygon_kwargs:
        polygon_kwargs["linewidth"] = 0.75 * matplotlib.rcParams["lines.linewidth"]
    if labelcolor is not None:
        label_kwargs["color"] = labelcolor
    if "color" not in label_kwargs:
        label_kwargs["color"] = polygon_kwargs["color"]
    if "fontsize" not in label_kwargs:
        label_kwargs["fontsize"] = 0.75 * matplotlib.rcParams["font.size"]

    if inverted:
        width_inches = -width_inches
    if zorder is None:
        zorder = 10

    # For more information on coordinate transformations in Matplotlib see
    # https://matplotlib.org/3.1.1/tutorials/advanced/transforms_tutorial.html

    # Convert the origin into figure coordinates in inches
    origin_disp = ax.transData.transform(origin)
    origin_dpi = fig.dpi_scale_trans.inverted().transform(origin_disp)

    # Obtain the bottom-right corner in data coordinates
    corner_dpi = origin_dpi + width_inches * numpy.array([1.0, 0.0])
    corner_disp = fig.dpi_scale_trans.transform(corner_dpi)
    corner = ax.transData.inverted().transform(corner_disp)

    (x1, y1) = (origin[0], origin[1])
    x2 = corner[0]

    # The width of the triangle in data coordinates
    width = x2 - x1
    # Compute offset of the slope
    log_offset = y1 / (x1 ** slope)

    y2 = log_offset * ((x1 + width) ** slope)
    height = y2 - y1

    # The vertices of the slope
    a = origin
    b = corner
    c = [x2, y2]

    # Draw the slope triangle
    X = numpy.array([a, b, c])
    triangle = plt.Polygon(X[:3,:], fill=False, zorder=zorder, **polygon_kwargs)
    ax.add_patch(triangle)

    # Convert vertices into display space
    a_disp = ax.transData.transform(a)
    b_disp = ax.transData.transform(b)
    c_disp = ax.transData.transform(c)

    # Figure out the center of the triangle sides in display space
    bottom_center_disp = a_disp + 0.5 * (b_disp - a_disp)
    bottom_center = ax.transData.inverted().transform(bottom_center_disp)

    right_center_disp = b_disp + 0.5 * (c_disp - b_disp)
    right_center = ax.transData.inverted().transform(right_center_disp)

    # Label alignment depending on inversion parameter
    va_xlabel = "top" if not inverted else "bottom"
    ha_ylabel = "left" if not inverted else "right"

    # Label offset depending on inversion parameter
    offset_xlabel = [0.0, -0.33 * label_kwargs["fontsize"]] if not inverted else [0.0, 0.33 * label_kwargs["fontsize"]]
    offset_ylabel = [0.33 * label_kwargs["fontsize"], 0.0] if not inverted else [-0.33 * label_kwargs["fontsize"], 0.0]

    # Draw the slope labels
    ax.annotate(r'$1$', bottom_center, xytext=offset_xlabel, textcoords='offset points', ha="center", va=va_xlabel, zorder=zorder, **label_kwargs)
    ax.annotate(fr'${slope}$', right_center, xytext=offset_ylabel, textcoords='offset points', ha=ha_ylabel, va="center", zorder=zorder, **label_kwargs)


# Kratos results (4,4)
eigvals_quad_OSS_4_4 = [ 1888.62, 1888.62, 2227.16, 3606.27, 6418.51, 6418.51, 11406.6, 11406.6, 12800, 13639.9, 13639.9, 15065, 15065.1, 22220.7, 22220.7, 23520.3, 25072.4, 25072.4, 34666.9, 34666.9 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_4_4 = [ 2227.16, 4034.68, 4034.68, 6532.24, 6532.24, 7659.96, 12800, 15548, 15548, 17741.2, 17741.2, 20861.1, 20861.1, 23017.2, 23017.2, 25935.5, 25935.5, 32137.1, 36827, 36827 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_irr_4_4 = [ 2189.71, 3115.99, 3115.99, 6119.61, 6358, 6358, 12000, 14400, 14400, 15085.7, 15085.7, 16996.7, 16996.7, 22316, 22316, 23813.4, 23813.4, 26400, 33600, 33600 ]

# Kratos results (8,8)
eigvals_quad_OSS_8_8 = [ 2034.12, 2266.17, 2266.17, 3897.89, 5288.7, 5288.7, 7573.37, 7573.37, 8908.66, 9039.43, 9039.43, 11225.4, 11225.4, 14444.9, 15443.7, 15443.7, 19618, 19618, 20686.7, 20686.7 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_8_8 = [ -6813.26, -6813.26, 2034.12, 5298.81, 5298.81, 8800.97, 8908.66, 11294.8, 11294.8, 15539.7, 15539.7, 16138.7, 16138.7, 19610.2, 19610.2, 20868.2, 20868.2, 23334.5, 26129, 26129 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_8_8 = eigvals_quad_ASGS_8_8[2:]
eigvals_quad_irr_8_8 = [ 2025.44, 2999.12, 2999.12, 5262.53, 5262.53, 5972.22, 8758.82, 11207.6, 11207.6, 12464, 12464, 15179.9, 15179.9, 15354.9, 15354.9, 20699.6, 20699.6, 22475.3, 24478.4, 25432 ]

# Kratos results (16,16)
eigvals_quad_OSS_16_16 = [ 1988.78, 2709.38, 2709.38, 5012.84, 5023.07, 5023.07, 8136.47, 9067.82, 9067.82, 10207.6, 10207.6, 10775.4, 10775.4, 13453.4, 13453.4, 15597.9, 17735.7, 17735.7, 17801.2, 17801.2 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_16_16 = [ -27253, -27253, 1988.78, 2773.7, 2773.7, 4796.82, 5023.46, 5023.46, 8136.47, 10211.3, 10211.3, 13458.8, 13458.8, 17752.1, 17752.1, 19010.2, 21195.2, 21195.2, 27080.3, 27080.3 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_16_16 = eigvals_quad_ASGS_16_16[2:]
eigvals_quad_irr_16_16 = [ 1986.65, 2970.41, 2970.41, 5014.74, 5014.74, 5934.43, 8101.77, 10191.4, 10191.4, 11996.5, 11996.5, 13379.2, 13379.2, 14941.1, 14941.1, 17715.9, 17715.9, 18828.8, 21050.1, 21050.1 ]

# Kratos results (32,32)
eigvals_quad_OSS_32_32 = [ 1977.62, 2890.54, 2890.54, 4956.78, 4956.78, 5643.24, 7955.12, 9953.89, 9953.89, 10837.7, 10837.7, 12984.4, 12984.4, 13256.5, 13256.5, 17016.9, 17016.9, 18067.8, 20052, 20092.3 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_32_32 = [ 1977.62, 2923.37, 2923.37, 4956.78, 4956.78, 5753.2, 7955.12, 9953.82, 9953.82, 11094.8, 11094.8, 12148.4, 12148.4, 12989.6, 12989.6, 13502.1, 13502.1, 17017.8, 17017.8, 18067.8 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_irr_32_32 = [ 1977.09, 2963.26, 2963.26, 4954.66, 4954.66, 5924.93, 7946.61, 9949.18, 9949.18, 11881.6, 11881.6, 12965.3, 12965.3, 14838.5, 14838.5, 17008.8, 17008.8, 18024.4, 20059, 20059 ]

# Kratos results (64,64)
eigvals_quad_OSS_64_64 = [ 1974.85, 2942.77, 2942.77, 4940.29, 4940.29, 5848.05, 7910.5, 9890.63, 9890.63, 11562.2, 11562.2, 12868.8, 12868.8, 14361.2, 14361.2, 16837.8, 16837.8, 17840.4, 19827.1, 19827.1 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_64_64 = [ 1974.85, 2951.95, 2951.95, 4940.29, 4940.29, 5883.5, 7910.5, 9890.63, 9890.63, 11693.5, 11693.5, 12868.8, 12868.8, 14558.1, 14558.1, 16837.8, 16837.8, 17840.4, 19827.1, 19827.1 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_irr_64_64 = [ 1974.71, 2961.48, 2961.48, 4939.76, 4939.76, 5922.56, 7908.38, 9889.44, 9889.44, 11853, 11853, 12864, 12864, 14812.9, 14812.9, 16835.7, 16835.7, 17829.6, 19818.6, 19818.6 ]

# Kratos results (128,128)
eigvals_quad_OSS_128_128 = [ 1974.15, 2956.32, 2956.32, 4936.17, 4936.17, 5903.06, 7899.38, 9874.86, 9874.86, 11771.1, 11771.1, 12840.1, 12840.1, 14689.5, 14689.5, 16793.2, 16793.2, 17784, 19761.2, 19761.2 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_128_128 = [ -10194.1, -10194.1, 1974.15, 2958.68, 2958.68, 3397.19, 3397.19, 4936.17, 4936.17, 5912.41, 7899.38, 9874.86, 9874.86, 11807.8, 11807.8, 12840.1, 12840.1, 14746.4, 14746.4, 16793.2 ] #tau(s) from paper c1=1 c2=4
eigvals_quad_ASGS_128_128 = eigvals_quad_ASGS_128_128[2:]
eigvals_quad_irr_128_128 = [ 1974.12, 2961.03, 2961.03, 4936.04, 4936.04, 5921.96, 7898.86, 9874.56, 9874.56, 11845.9, 11845.9, 12838.9, 12838.9, 14806.5, 14806.5, 16792.7, 16792.7, 17781.3, 19759, 19759 ]

# Overkill mesh (1250,1250) reference results (irreducible quad)
eigvals_ref = [ 1973.92, 2960.88, 2960.88, 4934.82, 4934.82, 5921.76, 7895.72, 9869.66, 9869.66, 11843.6, 11843.6, 12830.6, 12830.6, 14804.4, 14804.4, 16778.5, 16778.5, 17765.5, 19739.4, 19739.4 ]

# Calculate the errors with respect to the overkill mesh
eigvals_tri_OSS_err = []
eigvals_quad_OSS_err = []
eigvals_tri_ASGS_err = []
eigvals_quad_ASGS_err = []
eigvals_tri_irr_err = []
eigvals_quad_irr_err = []

eigval_id = 3
len_side = 1.0
n_div_vect = [4, 8, 16, 32, 64, 128]
h_vect = [len_side/n_div for n_div in n_div_vect]

for n_div in n_div_vect:
    # tri_OSS = locals()[f"eigvals_tri_OSS_{n_div}_{n_div}"][eigval_id]
    quad_OSS = locals()[f"eigvals_quad_OSS_{n_div}_{n_div}"][eigval_id]
    # tri_ASGS = locals()[f"eigvals_tri_ASGS_{n_div}_{n_div}"][eigval_id]
    quad_ASGS = locals()[f"eigvals_quad_ASGS_{n_div}_{n_div}"][eigval_id]
    # tri_irr = locals()[f"eigvals_tri_irr_{n_div}_{n_div}"][eigval_id]
    quad_irr = locals()[f"eigvals_quad_irr_{n_div}_{n_div}"][eigval_id]

    # eigvals_tri_OSS_err.append(abs(eigvals_ref[eigval_id] - tri_OSS))
    eigvals_quad_OSS_err.append(abs(eigvals_ref[eigval_id] - quad_OSS) / abs(eigvals_ref[eigval_id]))
    # eigvals_tri_ASGS_err.append(abs(eigvals_ref[eigval_id] - tri_ASGS))
    eigvals_quad_ASGS_err.append(abs(eigvals_ref[eigval_id] - quad_ASGS) / abs(eigvals_ref[eigval_id]))
    # eigvals_tri_irr_err.append(abs(eigvals_ref[eigval_id] - tri_irr))
    eigvals_quad_irr_err.append(abs(eigvals_ref[eigval_id] - quad_irr) / abs(eigvals_ref[eigval_id]))

## Calculate convergence Parameters
err_linear = []
err_quadratic = []
err_0_linear = 1
err_0_quadratic = 0.1
for h in h_vect:
    err_linear.append(err_0_linear*h)
    err_quadratic.append(err_0_quadratic*h**2)

## Plot the computed pressure coefficients solutions
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 3

locations_list = ['lower right','lower right','upper right','upper right', 'best']
linewidth = 1.0
legend_fontsize = 14
math_label_fontsize = 14

# Plot data
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(h_vect, eigvals_quad_OSS_err, '-s', color='k', linewidth=linewidth, markersize=4, alpha=0.9, label=r'OSS')
ax.loglog(h_vect, eigvals_quad_ASGS_err, '-s', color='tab:blue', linewidth=linewidth, markersize=4, alpha=0.9, label=r'ASGS')
ax.loglog(h_vect, eigvals_quad_irr_err, '-s', color='tab:orange', linewidth=linewidth, markersize=4, alpha=0.9, label=r'Irr.')
# ax.loglog(h_vect, eigvals_tri_OSS_err, '-^', color='k', linewidth=linewidth, markersize=4, alpha=0.9, label=r'OSS')
# ax.loglog(h_vect, eigvals_tri_ASGS_err, '-^', color='tab:blue', linewidth=linewidth, markersize=4, alpha=0.9, label=r'ASGS')
# ax.loglog(h_vect, eigvals_tri_irr_err, '-^', color='tab:orange', linewidth=linewidth, markersize=4, alpha=0.9, label=r'Irr.')
ax.loglog(h_vect, err_linear, '--', color='k', linewidth=linewidth, markersize=4, alpha=0.7, label='')
ax.loglog(h_vect, err_quadratic, '-.', color='k', linewidth=linewidth, markersize=4, alpha=0.7, label='')

# Set legend
plt.legend(numpoints=1, markerscale=2, loc='best',fontsize=legend_fontsize, fancybox=False, edgecolor='k')

# Set plot limits
# plt.xlim(-0.11, 0.11)
plt.ylim(1.0e-4, 1.0e-1)

# Set labels
x_label = r'$h$'
y_label = r'$\left\lvert \lambda_{0} - \bar{\lambda}_{0} \right\rvert$'
plt.xlabel(x_label, fontsize = math_label_fontsize)
plt.ylabel(y_label, fontsize = math_label_fontsize)

# Set axis ticks
# plt.xticks(h_vect, [rf'${round(h,2)}$' for h in h_vect], fontsize=math_label_fontsize)
plt.xticks(fontsize=math_label_fontsize)
plt.yticks(fontsize = math_label_fontsize)

# Set brackground grid
# plt.grid()

# Draw convergence triangles
tri_width = 0.35
tri_kwargs = {
    "color" : 'k',
    "alpha" : 0.7
}
tri_label_kwargs = {
    "color" : 'k',
    "alpha" : 0.7
}
tri_orig_linear = [0.05,250.0]
tri_orig_quadratic = [0.0375,0.35]
draw_loglog_slope(fig, ax, tri_orig_linear, tri_width, 1, inverted=True, polygon_kwargs=tri_kwargs, label_kwargs=tri_label_kwargs)
draw_loglog_slope(fig, ax, tri_orig_quadratic, tri_width, 2, inverted=False, polygon_kwargs=tri_kwargs, label_kwargs=tri_label_kwargs)

# Save plot and clear
plt.tight_layout()
plt.savefig(f"eigenvalue_convergence_compressible_slip_{eigval_id}.png",format='png', dpi=300)
plt.clf()
