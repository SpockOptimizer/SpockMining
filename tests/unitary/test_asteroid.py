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
Unitary test: Asteroid class

Functions contained:
    - TBC
"""

#######################################################################
# Imports area
#######################################################################

# Generic / Built-in


# Other Libs
import pykep as pk

# Own Libs
from src.asteroid import Asteroid
from src import constants


# First asteroid into data/candidates.txt file
data = [0.000000000000000000e+00,
       2.884773793646187973e+10,
       1.256922153715406233e-03,
       1.513348659383001449e-01,
       2.931641944496741203e+00,
       4.190976524071941434e+00,
       1.404363935586369294e+00,
       9.218574858087470458e-01,
       0.000000000000000000e+00]


def test_parse_asteroid_data():
    """
    Unitary test for parsing asteroid data.
    """

    # Expected planet object
    expectedPlanetObj = pk.planet.keplerian(
        pk.epoch_from_iso_string(constants.ISO_T_START),
        (
            data[1],
            data[2],
            data[3],
            data[4],
            data[5],
            data[6],
        ),
        constants.MU_TRAPPIST,
        constants.G * data[7],
        1,
        1.1,
        "Asteroid " + str(int(data[0])),
    )

    # Get Keplerian elements from incoming data
    expectedKeplerianElements = [
        data[1],
        data[2],
        data[3],
        data[4],
        data[5],
        data[6],
    ]



    # Pass data to Asteroid class
    asteroid = Asteroid(data)

    assert asteroid.asteroidId == data[0]
    assert asteroid.keplerianElements == expectedKeplerianElements
    assert asteroid.normalizedMass == data[7]
    assert asteroid.materialType == "Gold"

    # Test coordinates at t0
    expectedX = -18077450527.292923
    expectedY = 22304245037.175278
    expectedZ = -2752170398.901616

    x,y,z = asteroid.get_coordinates()
    assert x == expectedX
    assert y == expectedY
    assert z == expectedZ