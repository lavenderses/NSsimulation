import numpy as np
"""Plot Evry Particle.
"""


def plot_particle(ux, uy, ps, delt, delx, dely, ax):
    """Plot particles.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray): Air velocity in y dimension.
        ps (np.ndarray): Pressure.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        ax (matplotlib.pyplot.axes): Ax object prepared in main function.
    """
    p_map = np.copy(ps)
    p_map[:, 0] /= delx
    p_map[:, 1] /= dely
    p_map = p_map.astype(np.int)
    p_map[:, 0] = np.clip(p_map[:, 0], 0, uy.shape[0] - 1)
    p_map[:, 1] = np.clip(p_map[:, 1], 0, ux.shape[1] - 1)

    for i in range(len(ps)):
        p = p_map[i]
        ps[i, 0] += ux[p[0], p[1]] * delt
        ps[i, 1] += uy[p[0], p[1]] * delt

    ax.scatter(ps[:, 0], ps[:, 1], c='pink')

    return
