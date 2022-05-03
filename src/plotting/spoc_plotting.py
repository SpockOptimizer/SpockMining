# -*- coding: utf-8 -*-
#
#   Spock mining challenge
#   Copyright (C) 2022  Óscar Criado, Carlos Moreno
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#   See LICENSE

"""
Plotting module

Functions contained:
    - plotting_asteroids_by_material
    - plotting_zenith_by_material
    - plot_asteroids
    - plot_zenith_asteroids
"""
#######################################################################
# Imports area
#######################################################################

# Generic / Built-in


# Other Libs
import numpy as np
import matplotlib.pyplot as plt
import pykep as pk

# Own Libs
from src.asteroid import Asteroid


#######################################################################

def plotting_asteroids_by_material(asteroidList: list,
                                   materialType: str,
                                   ax: plt.Axes = None,
                                   color: str = 'k',
                                   alpha: float = 0.5,
                                   size: int = 5):
    """
    Plot asteroids by material type

    Parameters
    ----------
    asteroidList : list
        List of asteroids

    materialType : str
        Material type to plot

    ax : plt.Axes
        Axes to plot on

    color : str
        Color to plot

    alpha : float
        Alpha value to plot

    size : int
        Size to plot

    Returns
    -------
    None

    """

    asteroidsByMaterial = [asteroid for asteroid in asteroidList \
        if asteroid.materialType == materialType]

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    xCoord = []
    yCoord = []
    zCoord = []
    for asteroid in asteroidsByMaterial:
        x, y, z = asteroid.get_coordinates()
        xCoord.append(x)
        yCoord.append(y)
        zCoord.append(z)

    ax.scatter(xCoord, yCoord, zCoord,
               color=color,
               alpha=alpha,
               label=materialType,
               s=size)

def plotting_zenith_by_material(asteroidList: list,
                                materialType: str,
                                ax: plt.Axes = None,
                                color: str = 'k',
                                alpha: float = 0.5,
                                size: int = 5):
    """
    Zenith asteroids by material type

    Parameters
    ----------
    asteroidList : list
        List of asteroids

    materialType : str
        Material type to plot

    ax : plt.Axes
        Axes to plot on

    color : str
        Color to plot

    alpha : float
        Alpha value to plot

    size : int
        Size to plot

    Returns
    -------
    None

    """

    asteroidsByMaterial = [asteroid for asteroid in asteroidList \
        if asteroid.materialType == materialType]

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    xCoord = []
    yCoord = []
    for asteroid in asteroidsByMaterial:
        x, y, _ = asteroid.get_coordinates()
        xCoord.append(x)
        yCoord.append(y)

    ax.scatter(xCoord, yCoord,
               color=color,
               alpha=alpha,
               label=materialType,
               s=size)


def plot_asteroids(asteroidList: list,
                   figurename: str = 'asteroids.png'):
    """
    Plot all asteroids in list at the same plot

    Parameters
    ----------
    asteroidList : list
        List of asteroids

    figurename : str
        Name of the figure to save

    Returns
    -------
    None

    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plotting_asteroids_by_material(asteroidList,
                                   'Gold',
                                   ax,
                                   color='r')
    plotting_asteroids_by_material(asteroidList,
                                   'Platinum',
                                   ax,
                                   color='g')
    plotting_asteroids_by_material(asteroidList,
                                   'Nickel',
                                   ax,
                                   color='b')
    plotting_asteroids_by_material(asteroidList,
                                   'Propellant',
                                   ax,
                                   color='c')

    handles, labels = plt.gca().get_legend_handles_labels()
    byLabel = dict(zip(labels, handles))

    # Set title
    ax.title.set_text("Asteroids")

    # Set legend box below the plot
    ax.legend(byLabel.values(), byLabel.keys(),
                loc='lower left',
                bbox_to_anchor=(0, -0.1, 1, -0.1),
                ncol=5,
                mode="expand",
                borderaxespad=0.,
                )

    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Save figure
    plt.savefig(figurename)

def plot_zenith_asteroids(asteroidList: list,
                          figurename: str = 'zenith_asteroids.png'):
    """
    Plot zenith view of asteroids

    Parameters
    ----------
    asteroidList : list
        List of asteroids

    figurename : str
        Name of the figure to save

    Returns
    -------
    None

    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot zenith view by material type
    plotting_zenith_by_material(asteroidList,
                                'Gold',
                                ax,
                                color='r')
    plotting_zenith_by_material(asteroidList,
                                'Platinum',
                                ax,
                                color='g')
    plotting_zenith_by_material(asteroidList,
                                'Nickel',
                                ax,
                                color='b')
    plotting_zenith_by_material(asteroidList,
                                'Propellant',
                                ax,
                                color='c')

    # Avoid duplicated legend
    handles, labels = plt.gca().get_legend_handles_labels()
    byLabel = dict(zip(labels, handles))

    # Set aspect ratio to 1:1
    ax.set_aspect('equal')

    # Set grid
    ax.grid(True)
    ax.grid(which='minor', alpha=0.2)

    # Set title
    ax.title.set_text("Zenith view of asteroids")

    # Set legend box below the plot
    ax.legend(byLabel.values(), byLabel.keys(),
                loc='lower left',
                bbox_to_anchor=(0, -0.1, 1, -0.1),
                ncol=5,
                mode="expand",
                borderaxespad=0., )

    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # Save figure
    plt.savefig(figurename)

if __name__ == '__main__':

    # Load asteroids
    lines = np.loadtxt('data/candidates.txt')
    asteroids = [Asteroid(line) for line in lines]

    # Plot asteroids
    plot_asteroids(asteroids)
    plot_zenith_asteroids(asteroids)
