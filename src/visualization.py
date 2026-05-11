#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt  

mymark = ["o", "^", "x", "D", "+", "v", "<", ">", "h", "H"]

mycolor = [
    "#AD002D", "#1e50a2", "#69821b", "#f055f0", "#afafb0", 
    "#0095b9", "#89c3eb", "red", "blue", "orange", "pink"] 

def plotstyle(n):
    """
    Return n-th plot style.
    """
    n_col = len(mycolor)
    n_mark = len(mymark)
    n_total = n_col*n_mark
    assert n < n_total, f"Large index detected {n}."

    color = mycolor[int(n%n_col)]
    marker = mymark[int(np.floor(n/n_col))]
    return color, marker


def figure4lc(nline):
    """
    Create and return fig and axes (6x4) for lightcurve plot.
    If nline < 6, return a small sized figure.
    """
    # Absolute values
    N_raw            = 4
    N_line           = nline
    figsize_per_line = 2.5
    figsize_xlabel   = 0.05
    figsize_line     = figsize_per_line*N_line + figsize_xlabel
    figsize_raw      = 10
    print(f"N_raw, N_line = {N_raw}, {N_line}")

    # Relative values
    ## space for x label
    offset_xlabel = (figsize_xlabel/figsize_line)

    # Offset
    offset_width  = 0.0625
    offset_height = 0.03*(6/N_line)

    raw_width   = 1/N_raw*0.98
    line_height = 1/N_line*0.98

    # Figure size in each region
    fig_width  = raw_width*0.7
    fig_height = line_height*0.7

    axes = []
    n_fig = 0
    
    fig = plt.figure(figsize=(figsize_raw, figsize_line))
    for idx_line, line in enumerate(range(N_line)):
        for idx_raw, raw in enumerate(range(N_raw)):
            x0 = offset_width + idx_raw*raw_width
            y0 = offset_xlabel + offset_height + (N_line-idx_line-1)*line_height
            x1 = fig_width 
            y1 = fig_height
            ax = fig.add_axes([x0, y0, x1, y1]) 
            axes.append(ax)
            n_fig += 1
        if idx_line == (N_line-1):
            return fig, axes
    return fig, axes


def plot_lc_wmodel(df, JD0, ylim=None, label=None, out="lc.txt"):
    """Plot lightcurves in a 3x3 grid.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe.
    JD0 : float
        Zero Julian Day.
    ylim : tuple, optional
        y-axis limits.
    label : array-like, optional
        label
    out : str
        Output filename.
    """
    df = df.reset_index(drop=True)
    n_lc_list = list(set(df.n_lc))
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    axes = axes.flatten()

    for idx, ax in enumerate(axes):
        if idx >= len(n_lc_list):
            ax.axis("off")  # Hide unused subplot
            continue

        df_temp = df[df["n_lc"] == n_lc_list[idx]]
        ax.set_xlabel(f"JD-{JD0} [day]")
        ax.set_ylabel("Relative flux")
        if label:
            labelobs = f"Obs {label[idx]}"
        else:
            labelobs = f"Obs {idx+1}"
        labelmod = f"Model"
        ax.scatter(
            df_temp["jd"] - JD0, df_temp["flux"],
            color="black", label=labelobs, s=10)
        ax.scatter(
            df_temp["jd"] - JD0, df_temp["flux_model"],
            marker="x", color="red", label=labelmod)
        ax.legend()

    if ylim:
        ymin, ymax = ylim
        for ax in axes:
            ax.set_ylim([ymin, ymax])

    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def plot_plc_wmodel(df, JD0, rotP, Pishour, ylim=None, label=None, out="plc.png"):
    """Plot phase lightcurves in a 3x3 grid.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe
    JD0 : float
        Zero Julian day
    rotP : float
        Rotation period (in hours if Pishour=True, else in seconds)
    Pishour : bool
        Whether rotation period is in hours
    ylim : tuple
        (ymin, ymax) for y-axis range
    label : array-like, optional
        label
    out : str
        Output filename
    """
    import matplotlib.pyplot as plt

    df = df.reset_index(drop=True)
    n_lc_list = sorted(set(df.n_lc))  # sort for consistency

    # Convert rotation period
    rotP_sec = rotP * 3600. if Pishour else rotP
    rotP_day = rotP_sec / 86400.

    # Compute phase
    df["phase"] = ((df["jd"] - JD0) / rotP_day) % 1

    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    axs = axs.flatten()

    for idx, ax in enumerate(axs):
        if idx >= len(n_lc_list):
            ax.axis('off')  # Hide unused subplot
            continue
        df_temp = df[df["n_lc"] == n_lc_list[idx]]
        ax.set_xlabel("Rotational Phase")
        ax.set_ylabel("Relative flux")
        if label:
            labelobs = f"Obs {label[idx]}"
        else:
            labelobs = f"Obs {idx+1}"
        labelmod = f"Model"
        ax.scatter(
            df_temp["phase"], df_temp["flux"], color="black", s=10, label=labelobs)
        ax.scatter(
            df_temp["phase"], df_temp["flux_model"], marker="x", color="red", label=labelmod)
        ax.set_xlim([0.0, 1.0])
        ax.legend()

    if ylim:
        ymin, ymax = ylim
        for ax in axs:
            ax.set_ylim([ymin, ymax])

    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()



# Ellipsoid plot ==============================================================
# https://github.com/aleksandrbazhin/ellipsoid_fit_python
# https://jp.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
def fit_ellipsoid(x, y, z, align=False):
    # fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx +
    # 2Hy + 2Iz + J = 0 and A + B + C = 3 constraint removing one extra
    # parameter
    D = np.array([x * x + y * y - 2 * z * z,
                 x * x + z * z - 2 * y * y,
                 2 * x * y,
                 2 * x * z,
                 2 * y * z,
                 2 * x,
                 2 * y,
                 2 * z,
                 1 - 0 * x])
    d2 = np.array(x * x + y * y + z * z).T 
    u = np.linalg.solve(D.dot(D.T), D.dot(d2))
    a = np.array([u[0] + 1 * u[1] - 1])
    b = np.array([u[0] - 2 * u[1] - 1])
    c = np.array([u[1] - 2 * u[0] - 1])
    v = np.concatenate([a, b, c, u[2:]], axis=0).flatten()


    if align:
        # fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Gx + 2Hy + 2Iz = 1
        D = np.array(
            [x*x + y*y - 2*z*z,
            x*x + z*z - 2*y*y,
              2*x,
              2*y,
              2*z,
            1 + 0*x])
        d2 = np.array(x*x + y*y + z*z).T 
        u = np.linalg.solve(D.dot(D.T), D.dot(d2))
        a = np.array([u[0] + 1 * u[1] - 1])
        b = np.array([u[0] - 2 * u[1] - 1])
        c = np.array([u[1] - 2 * u[0] - 1])
        d = np.array([0, 0, 0])
        v = np.concatenate([a, b, c, d, u[2:6]], axis=0).flatten()


    A = np.array([[v[0], v[3], v[4], v[6]],
                  [v[3], v[1], v[5], v[7]],
                  [v[4], v[5], v[2], v[8]],
                  [v[6], v[7], v[8], v[9]]])

    center = np.linalg.solve(- A[:3, :3], v[6:9])

    translation_matrix = np.eye(4)
    translation_matrix[3, :3] = center.T

    R = translation_matrix.dot(A).dot(translation_matrix.T)

    evals, evecs = np.linalg.eig(R[:3, :3] / -R[3, 3])
    evecs = evecs.T

    radii = np.sqrt(1. / np.abs(evals))
    radii *= np.sign(evals)

    return center, evecs, radii, v


# https://github.com/minillinim/ellipsoid
def draw_ellipsoid(center, radii, rotation, ax, plot_axes=False, cage_color='b', cage_alpha=0.2):
    """Plot an ellipsoid"""
        
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    
    # cartesian coordinates that correspond to the spherical angles:
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    # rotate accordingly
    for i in range(len(x)):
        for j in range(len(x)):
            [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], rotation) + center

    if plot_axes:
        # make some purdy axes
        axes = np.array([[radii[0],0.0,0.0],
                         [0.0,radii[1],0.0],
                         [0.0,0.0,radii[2]]])
        # rotate accordingly
        for i in range(len(axes)):
            axes[i] = np.dot(axes[i], rotation)

        # plot axes
        for p in axes:
            X3 = np.linspace(-p[0], p[0], 100) + center[0]
            Y3 = np.linspace(-p[1], p[1], 100) + center[1]
            Z3 = np.linspace(-p[2], p[2], 100) + center[2]
            ax.plot(X3, Y3, Z3, color=cage_color)

    # plot ellipsoid
    ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=cage_color, alpha=cage_alpha)

# 
# 
# def set_axes_equal(ax: plt.Axes):
#     ax.set_box_aspect([1,1,1])
#     limits = np.array([
#         ax.get_xlim3d(),
#         ax.get_ylim3d(),
#         ax.get_zlim3d(),
#     ])
#     x, y, z = np.mean(limits, axis=1)
#     radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
#     ax.set_xlim3d([x - radius, x + radius])
#     ax.set_ylim3d([y - radius, y + radius])
#     ax.set_zlim3d([z - radius, z + radius])
# 
# 
# if __name__=='__main__':
#     data = np.loadtxt("mag_out.txt")
#     data2 = data_regularize(data, divs=8)
# 
#     center, evecs, radii, v = ellipsoid_fit(data2)
# 
#     data_centered = data - center.T
#     data_centered_regularized = data2 - center.T
# 
#     a, b, c = radii
#     r = (a * b * c) ** (1. / 3.)
#     D = np.array([[r/a, 0., 0.], [0., r/b, 0.], [0., 0., r/c]])
#     #http://www.cs.brandeis.edu/~cs155/Lecture_07_6.pdf
#     #affine transformation from ellipsoid to sphere (translation excluded)
#     TR = evecs.dot(D).dot(evecs.T)
#     data_on_sphere = TR.dot(data_centered_regularized.T).T
# 
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
# 
#     # for direction in (-1, 1):
#     #     for point in np.diag(direction * np.max(data) * np.array([1, 1, 1])):
#     #         ax.plot([point[0]], [point[1]], [point[2]], 'w')
# 
#     ax.scatter(data_centered[:,0], data_centered[:,1], data_centered[:,2], marker='o', color='g')
#     # ax.scatter(data_centered_regularized[:, 0], data_centered_regularized[:, 1],
#     #            data_centered_regularized[:, 2], marker='o', color='b')
#     ax.scatter(data_on_sphere[:, 0], data_on_sphere[:, 1],
#                data_on_sphere[:, 2], marker='o', color='r')
# 
#     ellipsoid_plot([0, 0, 0], [r, r, r], evecs, ax=ax, plot_axes=True, cage_color='orange')
# 
#     #ax.plot([r],[0],[0],color='r',marker='o')
#     #ax.plot([radii[0]],[0],[0],color='b',marker='o')
# 
#     set_axes_equal(ax)
#     plt.show()
# 
# Ellipsoid plot ==============================================================
