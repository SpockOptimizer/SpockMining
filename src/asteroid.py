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
Asteroid class definition file
"""
#######################################################################
# Imports area
#######################################################################

# Generic / Built-in


# Other Libs
import pykep as pk

# Own Libs
from src import constants


#######################################################################


class Asteroid:
    """
    Asteroid class

    Attributes
    ----------
    asteroidDaa : list
        Asteroid data

    Methods
    -------
    parse_asteroid_data
    get_coordinates

    """

    def __init__(self,
                 asteroidData: str):
        """
        Initialize asteroid object
        """
        self.parse_asteroid_data(asteroidData)

    def parse_asteroid_data(self,
                            asteroidData: list):
        """
        Parse candidate asteroid data

        Parameters
        ----------
        asteroidData : list
            Asteroid data

        """
        self.asteroidId = asteroidData[0]
        self.keplerianElements = [asteroidData[1],
                                   asteroidData[2],
                                   asteroidData[3],
                                   asteroidData[4],
                                   asteroidData[5],
                                   asteroidData[6]]
        self.planetObject = pk.planet.keplerian(
            pk.epoch_from_iso_string(constants.ISO_T_START),
            (
                asteroidData[1],
                asteroidData[2],
                asteroidData[3],
                asteroidData[4],
                asteroidData[5],
                asteroidData[6],
            ),
            constants.MU_TRAPPIST,
            constants.G * asteroidData[7],
            1,
            1.1,
            "Asteroid " + str(int(asteroidData[0])),
        )

        # And asteroids' masses and material type
        self.normalizedMass = asteroidData[-2]
        self.materialType = set_material_type(int(asteroidData[-1]))

    def get_coordinates(self):
        """
        Get coordinates of asteroid
        """
        r, _ = self.planetObject.eph(
            pk.epoch_from_iso_string(constants.ISO_T_START))
        x, y, z = r
        return x, y, z

def set_material_type(materialType: int)->str:
    """
    Transform material type to string

    Parameters
    ----------
    materialType : int
        Material type integer to be transformed

    Returns
    -------
    str
        Material type string

    """
    if materialType == 0:
        return "Gold"
    if materialType == 1:
        return "Platinum"
    if materialType == 2:
        return "Nickel"
    return "Propellant"
