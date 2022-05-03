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
Module to store all constants used in the program
"""

# Start and end epochs
ISO_T_START = "30190302T000000"
ISO_T_END = "30240302T000000"

# Cavendish constant (m^3/s^2/kg)
G = 6.67430e-11

# Sun_mass (kg)
SM = 1.989e30

# Mass and Mu of the Trappist-1 star
MS = 8.98266512e-2 * SM
MU_TRAPPIST = G * MS

# DV per propellant [m/s]
DV_PER_PROPELLANT = 10000

# Maximum time to fully mine an asteroid
TIME_TO_MINE_FULLY = 30

# Eccentric anomaly constants
MAX_ITER_E = 15
EPSILON_E = 1e-10
