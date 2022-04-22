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
Conversor module

Functions contained:
    - coe_to_rv
    - mean_to_eccentric
"""
#######################################################################
# Imports area
#######################################################################

# Generic / Built-in


# Other Libs
import numpy as np

# Own Libs
import constants


#######################################################################


def coe_to_rv(μ: float, coe: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Calculate cartesian coordinates from keplerian orbital elements
    given a time

    Parameters
    ----------
    μ : float [m**3/s**2]
        Standard gravitational parameter

    coe : np.ndarray (6xN) [m, -, rad, rad, rad, rad]
        Matrix of keplerian elements [a, e, i, Ω, AoP, TA]

    t : np.ndarray (Mx1) [s]
        Array of times to compute de cartesian coordinates


    Returns
    -------
    positions : np.ndarray (NxMx3) [m]
        Cartesian positions for times t

    """
    sma, ecc, inc, Ω, omega, trueAnomaly = coe.T

    # Vector lenghts
    n = len(sma)
    m = len(t)
    transMat = np.ones((m, n))
    positions = np.zeros((n, m, 3))

    # Transform the coe to be able to operate with in form of matrices
    sma = sma * transMat
    ecc = ecc * transMat
    inc = inc * transMat
    Ω = Ω * transMat
    omega = omega * transMat
    trueAnomaly = trueAnomaly * transMat

    t = np.array([t]).T * transMat

    # Calculate true anomaly through time
    meanMotion = np.sqrt(μ/sma**3)
    trueAnomaly = trueAnomaly + meanMotion * t

    eccAnomaly = np.arctan2(np.sqrt(1 - ecc**2) * np.sin(trueAnomaly),
                            ecc + np.cos(trueAnomaly))

    # Calculate satellite angular position in the orbital plane
    Φ = trueAnomaly + omega

    # Calculate radiovector in the orbital plane
    r = sma*(1 - ecc*np.cos(eccAnomaly))

    # Calculate positions in the orbital plane
    xOrbPlane = r * np.cos(Φ)
    yOrbPlane = r * np.sin(Φ)

    # Calculate position in ECI
    positions[:, :, 0] = (xOrbPlane*np.cos(Ω)
                          - yOrbPlane*np.cos(inc)*np.sin(Ω)).T
    positions[:, :, 1] = (xOrbPlane*np.sin(Ω)
                          + yOrbPlane*np.cos(inc)*np.cos(Ω)).T
    positions[:, :, 2] = (yOrbPlane*np.sin(inc)).T

    return positions


def mean_to_eccentric(ecc: np.ndarray, meanAnomaly: np.ndarray) -> np.ndarray:
    """
    Calculate eccentric anomaly through an iterative method

    Parameters
    ----------
    ecc : float [-]
        Eccentricity

    meanAnomaly : np.ndarray (N) [rad]
        Mean anomaly

    Returns
    -------
    eccAnomaly : np.ndarray (N) [rad]
        Eccentric anomaly

    """
    eccAnomaly = meanAnomaly
    count = 0

    while count < constants.MAX_ITER_E:
        oldEccAnomaly = eccAnomaly

        eccAnomaly = \
            eccAnomaly - (eccAnomaly - ecc*np.sin(eccAnomaly) - meanAnomaly)\
            / (1 - ecc*np.cos(eccAnomaly))

        if count > constants.MAX_ITER_E:
            return eccAnomaly

        if np.all((eccAnomaly - oldEccAnomaly) <= constants.EPSILON_E):
            return eccAnomaly

        count += 1
