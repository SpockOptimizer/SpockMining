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
Unitary test: Reference system coordinates library

Functions contained:
    - coe2rv
    - calculate_eccentric_anomaly
"""
#######################################################################
# Imports area
#######################################################################

# Generic / Built-in


# Other Libs
import numpy as np

# Own Libs
from src.refsys import coordinates

#######################################################################


#######################################################################
# Constants area
#######################################################################

# e, M[deg], E[deg], nu[deg]
ELLIPTIC_ANGLES_DATA = [
    (0.0, 0.0, 0.0, 0.0),
    (0.01, 10.0, 10.10048, 10.20146),
    (0.05, 25.0, 26.26786, 27.56558),
    (0.10, 135.0, 138.77583, 142.42239),
    (0.20, 55.0, 65.42081, 76.37621)
    ]


############################
# MATH CONSTANTS
############################
PI = 3.1415926535898
DEG_TO_RAD = PI/180
RAD_TO_DEG = 180/PI


#######################################################################


def test_coe_to_rv():
    """
    Unitary test for calculating Keplerian elements to cartesian.

    References
    ----------
    Vallado, D.A., 2013. Fundamentals of Astrodynamicsand Applications,
    4th edition. Example 2-6, page 119
    """
    μ = 3.986004418e14
    coe = np.array([[36126642.834,
                     0.83285,
                     1.53362,
                     3.97743,
                     0.93166,
                     1.61155]])
    t = np.array([0])

    positions = np.array([[6525344, 6861535, 6449125]])

    assert np.allclose(coordinates.coe_to_rv(μ, coe, t), positions)


def test_mean_to_eccentric():
    """
    Unitary test for calculating eccentric anomaly from mean anomaly.

    References
    ----------
    Kent, J.T., Taack, G.B.,  Larson, D.C. & NASA, 1963. Tables
    for Eccentric and True Anomaly in Elliptic Orbits
    """

    for ecc, meanAnomaly, eccAnomaly, _ in ELLIPTIC_ANGLES_DATA:
        assert np.isclose(
            coordinates.mean_to_eccentric(ecc, meanAnomaly*DEG_TO_RAD)\
                *RAD_TO_DEG,
            eccAnomaly)
