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