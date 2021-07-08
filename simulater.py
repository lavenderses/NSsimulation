import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as anime
from decimal import Decimal, ROUND_HALF_UP
import functions as F2d
import functions3d as F3d
"""The Simulation of Navier Stokes Equatin in 2D & 3D.
"""


class Simulate2D:
    """2 dimensional simulated room model.

    Attributes:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        p (np.ndarray): Pessure.
        lx (float): Room length in x dimenstion.
        ly (float): Room length in y dimenstion.
        DELT (float): Infinitensimal time. Defaults to 0.01.
        DELL (float): Infinitensimal displacement. Defaults to 0.1.
        MU (float): Viscosity value. Defaults to 1.82e-5.
        EPS (float): Infinitensimal value for Jaccobi function.
                     Defaults to 1e-8.
        imgs (list): Axes list to animate.
        xx (np.ndarray): Cordinate array.
        yy (np.ndarray): Cordinate array.
        fig (matplotlib.pyplot.figure): Plot figure.
    """

    def __init__(self, lx, ly, DELT=0.01, DELL=0.1, RHO=1.293, MU=1.82e-5, EPS=1e-8):
        """Initialize this class and set attributes.

        Args:
            lx (float): Room length in x dimenstion.
            ly (float): Room length in y dimenstion.
            DELT (float, optional): Infinitensimal time. Defaults to 0.01.
            DELL (float, optional): Infinitensimal displacement. Defaults to 0.1.
            RHO (float, optional): Density field. Defaults to 1.293.
            MU (float, optional): Viscosity value. Defaults to 1.82e-5.
            EPS (float, optional): Infinitensimal value for Jaccobi function.
                                   Defaults to 1e-8.
        """
        divlx = int(Decimal(str(lx / DELL)).quantize(Decimal('0'),
                                                          rounding=ROUND_HALF_UP))
        divly = int(Decimal(str(ly / DELL)).quantize(Decimal('0'),
                                                          rounding=ROUND_HALF_UP))
        self.ux = np.zeros((divlx + 3, divly + 2))
        self.uy = np.zeros((divlx + 2, divly + 3))
        self.p = np.zeros((divlx + 2, divly + 2))
        self.lx, self.ly = lx, ly
        self.DELT, self.DELL = DELT, DELL
        self.MU, self.RHO, self.EPS = MU, RHO, EPS
        self.imgs = []

        # Points
        self.xx = np.array([np.full(divly + 2, i) for i in range(divlx + 2)])
        self.yy = np.tile(np.arange(divly + 2), divlx + 2)

        plt.rcParams['font.size'] = 26
        self.fig = plt.figure(figsize=(20, 14), dpi=100)

    def update(self, hs, v0, theta, time_range, w=None, rad_range=None, CNT_MAX=10000, save=False, arrow=True):
        """Simulate this model.

        Args:
            hs (list): Fan's cordinate.
            v0 (float): Fan's initial velocity.
            theta (float): Fan's elevation angle.
            time_range (float): Simulation time.
            w (float, optional): Fan's ngular velocity. Defaults to None.
            rad_range (list, optional): Fan's angle range. Defaults to None.
            CNT_MAX (int, optional): Iteration value. Defaults to 10000.
            save (bool, optional): Whether you save every plot. Defaults to False.
            arrow (bool, optional): Plot style. Arrow or particle. Defaults to True.
        """
        lx, ly  = self.lx, self.ly
        DELT, DELL = self.DELT, self.DELL
        passed_time, num = 0, 0

        dirname = './imgs_ex/2d_{}_{}'.format(v0, theta)
        os.makedirs(dirname, exist_ok=True)
        theta = np.deg2rad(theta)

        # Particles Array
        ptclsx = np.random.rand(4000) * lx
        ptclsy = np.random.rand(4000) * ly
        particles = np.stack([ptclsx, ptclsy]).T

        while passed_time < time_range:
            self.ux, self.uy, theta, w = F2d.ff(self.ux, self.uy, self.DELT, hs, v0, theta, w, rad_range)

            # Canvus
            ax = self.fig.add_subplot(111,
                                      xlim=(-0.1, self.lx * 1.1),
                                      ylim=(-0.1, self.ly * 1.1),
                                      aspect='equal',)
            ax.set(xlabel='Room Length x/m', ylabel='Room Height y/m')
            ax.set_title('t = {:.2f} [s]'.format(num * DELT))

            # Plot
            if arrow:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes('right', '5%', pad='3%')
                im = F2d.plot(self.ux, self.uy, DELL, DELL, self.xx, self.yy, ax)
                self.fig.colorbar(im, cax=cax)
            else:
                im = F2d.pp(self.ux, self.uy, particles, DELT, DELL, DELL, ax)

            self.fig.savefig('{}/{:0=10}.png'.format(dirname, num))
            if not save:
                self.imgs.append([im])
            plt.cla()
            plt.clf()

            print('time : {:.2f}'.format(num * DELT))
            # Advection
            ux_ast = F2d.adv_x(self.ux, self.uy, DELT, DELL, DELL)
            uy_ast = F2d.adv_y(self.ux, self.uy, DELT, DELL, DELL)

            # Viscosity
            ux_ast, uy_ast = F2d.vis(self.ux, self.uy, ux_ast, uy_ast, DELL, DELL, self.MU, self.RHO)
    
            # Solce Poisson Eq to Align ux&uy (div u = 0)
            self.p = F2d.poisson(ux_ast, uy_ast, DELT, DELL, DELL, self.p, self.EPS, CNT_MAX)

            # Fix each Velocity
            self.ux, self.uy = F2d.fix_u(ux_ast, uy_ast, self.p, DELT, DELL, DELL, self.RHO)
    
            passed_time = num * DELT
            num += 1

        if not save:
            anim = anime.ArtistAnimation(self.fig, self.imgs, interval=DELT * 1000)
            anim.save('./animation_{}_{}.mp4'.format(v0, theta), writer='ffmpeg')
        

class Simulate3D:
    """3 dimensional simulated room model.

    Attributes:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        uz (np.ndarray)): Air velocity in z dimension.
        p (np.ndarray): Pessure.
        lx (float): Room length in x dimenstion.
        ly (float): Room length in y dimenstion.
        lz (float): Room length in z dimenstion.
        DELT (float): Infinitensimal time. Defaults to 0.01.
        DELL (float): Infinitensimal displacement. Defaults to 0.1.
        MU (float): Viscosity value. Defaults to 1.82e-5.
        RHO (float): Density field.
        EPS (float): Infinitensimal value for Jaccobi function.
                     Defaults to 1e-8.
        imgs (list): Axes list to animate.
        xx (np.ndarray): Cordinate array.
        yy (np.ndarray): Cordinate array.
        zz (np.ndarray): Cordinate array.
        fig (matplotlib.pyplot.figure): Plot figure.
    """

    def __init__(self, lx, ly, lz, DELT=0.01, DELL=0.1, RHO=1.293, MU=1.82e-5, EPS=1e-8):
        """Initialize this class and set attributes.

        Args:
            lx (float): Room length in x dimenstion.
            ly (float): Room length in y dimenstion.
            lz (float): Room length in z dimenstion.
            DELT (float, optional): Infinitensimal time. Defaults to 0.01.
            DELL (float, optional): Infinitensimal displacement. Defaults to 0.1.
            RHO (float, optional): Density field. Defaults to 1.293.
            MU (float, optional): Viscosity value. Defaults to 1.82e-5.
            EPS (float, optional): Infinitensimal value for Jaccobi function.
                                   Defaults to 1e-8.
        """
        # Set Constance Values
        divlx = int(Decimal(str(lx / DELL)).quantize(Decimal('0'),
                                                     rounding=ROUND_HALF_UP))
        divly = int(Decimal(str(ly / DELL)).quantize(Decimal('0'),
                                                     rounding=ROUND_HALF_UP))
        divlz = int(Decimal(str(lz / DELL)).quantize(Decimal('0'),
                                                     rounding=ROUND_HALF_UP))
        self.ux = np.zeros((divlx + 3, divly + 2, divlz + 2))
        self.uy = np.zeros((divlx + 2, divly + 3, divlz + 2))
        self.uz = np.zeros((divlx + 2, divly + 2, divlz + 3))
        self.p = np.zeros((divlx + 2, divly + 2, divlz + 2))
        self.lx, self.ly, self.lz = lx, ly, lz
        self.DELT, self.DELL = DELT, DELL
        self.MU, self.RHO, self.EPS = MU, RHO, EPS
        self.imgs = []
    
        self.xx = np.array([np.full((divlz + 2) * (divly + 2), i) for i in range(divlx + 2)])
        self.yy = np.tile(np.array([np.full(divlz + 2, j) for j in range(divly + 2)]).flatten(), divlx + 2)
        self.zz = np.tile(np.arange(divlz + 2), (divly + 2) * (divlx + 2))

        plt.rcParams['font.size'] = 26
        self.fig = plt.figure(figsize=(24, 14), dpi=100)

    def update(self, hs, v0, theta, phai, time_range, w=None, rad_range=None, CNT_MAX=10000, save=False, arrow=True):
        """Simulate this model.

        Args:
            hs (list): Fan's cordinate.
            v0 (float): Fan's initial velocity.
            theta (float): Fan's elevation angle.
            phai (float): Fan's azimuth angle.
            time_range (float): Simulation time.
            w (float, optional): Fan's ngular velocity. Defaults to None.
            rad_range (list, optional): Fan's angle range. Defaults to None.
            CNT_MAX (int, optional): Iteration value. Defaults to 10000.
            save (bool, optional): Whether you save every plot. Defaults to False.
            arrow (bool, optional): Plot style. Arrow or particle. Defaults to True.
        """
        lx, ly, lz  = self.lx, self.ly, self.lz
        DELT, DELL = self.DELT, self.DELL
        MU, RHO, EPS = self.MU, self.RHO, self.EPS
        passed_time, num = 0, 0

        dirname = './imgs_ex/3d_{}_{}_{}'.format(v0, theta, phai)
        os.makedirs(dirname, exist_ok=True)

        theta = np.deg2rad(theta)
        phai = np.deg2rad(phai)
        phai = np.deg2rad(phai)

        # Particles Array
        ptclsx = np.random.rand(4000) * lx
        ptclsy = np.random.rand(4000) * ly
        ptclsz = np.random.rand(4000) * lz
        particles = np.stack([ptclsx, ptclsy, ptclsz]).T

        while passed_time < time_range:
            # First Condition
            self.ux, self.uy, self.uz = F3d.ff(self.ux, self.uy, self.uz, DELT, hs, v0, theta, phai)

            # Canvus
            ax = self.fig.add_subplot(111,
                                      projection='3d',
                                      xlim=(-0.1, lx * 1.1),
                                      ylim=(-0.1, ly * 1.1),
                                      zlim=(-0.1, lz * 1.1),)
            ax.set_xlabel('Room Length x/m', labelpad=20)
            ax.set_ylabel('Room Depth y/m', labelpad=20)
            ax.set_zlabel('Room Height z/m', labelpad=20)
            ax.set_title('t = {:.2f} [s]'.format(num * DELT))

            # Plot
            if arrow:
                img = F3d.plot(self.ux, self.uy, self.uz, DELL, DELL, DELL, self.xx, self.yy, self.zz, ax)
                self.fig.colorbar(img)
            else:
                img = F3d.p_plot(self.ux, self.uy, self.uz, lx, ly, lz, particles, DELT, DELL, DELL, DELL, ax)

            self.fig.savefig('{}/{:0=10}.png'.format(dirname, num))
            if not save:
                self.imgs.append([img])
            plt.cla()
            plt.clf()

            print('time : {:.2f}'.format(num * DELT))
            # Advection
            ux_ast = F3d.adv_x(self.ux, self.uy, self.uz, DELT, DELL, DELL, DELL)
            uy_ast = F3d.adv_y(self.ux, self.uy, self.uz, DELT, DELL, DELL, DELL)
            uz_ast = F3d.adv_z(self.ux, self.uy, self.uz, DELT, DELL, DELL, DELL)

            # Viscosity
            ux_ast, uy_ast, uz_ast = F3d.vis(self.ux, self.uy, self.uz, ux_ast, uy_ast, uz_ast, DELL, DELL, DELL, MU, RHO)

            # Solce Poisson Eq to Align ux&uy (div u = 0)
            self.p = F3d.poisson(ux_ast, uy_ast, uz_ast, DELT, DELL, DELL, DELL, self.p, EPS, CNT_MAX)

            # Fix each Velocity
            self.ux, self.uy, self.uz = F3d.fix_u(ux_ast, uy_ast, uz_ast, self.p, DELT, DELL, DELL, DELL, RHO)

            passed_time = num * DELT
            num += 1

        if not save:
            anim = anime.ArtistAnimation(self.fig, self.imgs, interval=100)
            anim.save('./animation_{}_{}_{}.mp4'.format(v0, theta, phai), writer='ffmpeg')


if __name__ == '__main__':
    # Set Each Constance
    room_x = 4.0
    room_y = 2.0
    room_z = 2.5
    h = [[1, 13],
         [1, 14]]
    h_3d = [[1, 3, 3],
            [1, 3, 4],
            [1, 4, 3],
            [1, 4, 4]]
    time_range = 30.

    rng = np.deg2rad(60)
    rad_range = [-rng, rng]
    T = 18.0
    w = 2 * np.pi / T

    '''
    the = 0
    v0 = 5.0
    model = Simulate2D(room_x, room_y)
    model.update(h, v0, the, time_range, arrow=False)
    '''
    the = 60
    pha = 30
    v0 = 15.0
    time_range = 120.
    model = Simulate3D(room_x, room_y, room_z)
    model.update(h_3d, v0, the, pha, time_range, arrow=False)
    #'''
