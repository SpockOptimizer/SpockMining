''' Base source code of SpOC Mining ESA Challenge
    from: https://optimize.esa.int/challenge/spoc-mining/p/mine-the-belt'''
import numpy as np
import pykep as pk
import matplotlib.pyplot as plt

################
### Constants
################

# Start and end epochs
T_START = pk.epoch_from_iso_string("30190302T000000")
T_END = pk.epoch_from_iso_string("30240302T000000")

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

# Loading the asteroid data
data = np.loadtxt("data/candidates.txt")
asteroids = []
for line in data:
    p = pk.planet.keplerian(
        T_START,
        (
            line[1],
            line[2],
            line[3],
            line[4],
            line[5],
            line[6],
        ),
        MU_TRAPPIST,
        G * line[7],  # these variable are not relevant for this problem
        1,  # these variable are not relevant for this problem
        1.1,  # these variable are not relevant for this problem
        "Asteroid " + str(int(line[0])),
    )
    asteroids.append(p)

# And asteroids' masses and material type
asteroidMasses = data[:, -2]
asteroidMaterials = data[:, -1].astype(int)


def convert_to_chromosome(solution, checkCorrectness=True):
    """Convert a solution to a fixed-length chromosome to be
       compatible with the optimize and pygmo framework.

    Args:
        solution (list or np.array): The solution to be converted. Has to
            have format [t_arrival_0, ..., t_arrival_n, t_mining_0, ...,
            t_mining_n, ast_id_0, ..., ast_id_n]
        checkCorrectness (bool): If True, the function will check if the
            solution fulfills some validity checks (unique asteroids,
            solution length).

    Returns:
        np.array: The chromosome as required by optimize / pygmo.
    """
    N = len(solution) // 3  # number of asteroids in chromosome
    ast_ids = solution[2 * N :]  # ids of the visited asteroids

    # Check if the solution is valid
    if checkCorrectness:
        assert (
            len(solution) % 3 == 0
        ), '''Solution must be a multiple of 3 containing asteroid id,
            arrival time and time to mine. for each visit.'''

        assert (
            len(set(ast_ids)) - len(ast_ids) == 0
        ), "Asteroid IDs must be unique, can only visit each asteroid once."

    # The final chromosome will need to contain all asteroids, so we need to
    # add the asteroids to the chromosome that are not in the solution
    chromosome = np.zeros(30000, dtype=np.float64)

    # Set placeholder values for mining times and arrival times for irrelevant
    # chromosome entries
    chromosome[N:10000] = 0
    chromosome[10000 + N : 20000] = 0

    # Add time of arrivals and mining times
    chromosome[:N] = solution[:N]
    chromosome[10000 : 10000 + N] = solution[N : 2 * N]

    # Add the asteroids that are in the solution
    chromosome[20000 : 20000 + N] = ast_ids

    # Add the asteroids that are not in the solution.
    # There is the potential of a very rare edgecase where by conincidence
    # the next asteroid added this way could still be visited but this is
    # excessively unlikely
    ast_not_in_solution =\
        set(np.arange(10000)).symmetric_difference(set(ast_ids))
    chromosome[20000 + N :] = np.array(list(ast_not_in_solution))

    return chromosome


class BeltMiningUdp:
    """
    pygmo User Defined Problem (UDP) describing the optimization problem.
    https://esa.github.io/pygmo2/tutorials/coding_udp_simple.html explains
    what more details on UDPs.
    """

    def __init__(self, missionWindow):
        """Initialize the UDP.

        Args:
            missionWindow (list [float, float]): Bounds on the overall
            mission in days.
        """
        self.asteroids = asteroids
        self.asteroidMasses = asteroidMasses
        self.asteroidMaterialsTypes = asteroidMaterials
        self.missionWindow = missionWindow
        self.n = len(self.asteroids)
        self.MU = MU_TRAPPIST

    def get_bounds(self):
        """Get bounds for the decision variables.

        Returns:
            Tuple of lists: bounds for the decision variables.
        """
        lb = [self.missionWindow[0]] * self.n + [0] * self.n + [0] * self.n
        ub = [self.missionWindow[1]] * self.n + [60] * self.n + [self.n - 1] *\
             self.n
        return (lb, ub)

    def get_nix(self):
        """Get number of integer variables.

        Returns:
            int: number of integer variables.
        """
        return self.n

    def get_nic(self):
        """Get number of inequality constraints.

        Returns:
            int: number of inequality constraints.
        """
        # Inequality constraints are only set to all visiting epochs (except
        # the first)
        return self.n - 1

    def get_nec(self):
        """Get number of equality constraints.

        Returns:
            int: number of equality constraints.
        """
        # The only equality constraint is that each asteroid must be in the
        # list exactly once
        return 1

    def _evaluate_journey(self, chromose, verbose=False):
        # propellant level of the ship, cannot go below 0 or we abort
        propellant = 1
        # viable number of visited asteroids, will be computed
        visited = 0
        # number of asteroids in chromosome
        n = len(chromose) // 3
        # time at arrival of each asteroid in days
        time_at_arrival = chromose[:n]
        # how many days spent mining each asteroid
        time_spent_preparing = chromose[n : 2 * n]
        # last is propellant and will be disregarded for score
        material_prepared = [
            0
        ] * 4
        # list of lists of materials prepared at each time step (for plotting)
        material_prepared_at_t = (
            []
        )
        # ids of the visited asteroids
        ast_ids = chromose[2 * n :]
        if verbose:
            print('''ID\tt0\tPropellant \tDV \t  Material ID\t Prepared\t
            \tScore''')

        propellant_at_time = []

        # Lets compute the fitness
        for i in range(1, n):

            # Get indices of currently visited asteroid
            # and the previous one
            current_ast_id = int(ast_ids[i])
            previous_ast_id = int(ast_ids[i - 1])

            ###################### Step 1 #######################
            # Validate the transfer from asteroid i to i+1      #
            #####################################################

            # Break as soon as we exceed mission window
            if time_at_arrival[i] - time_at_arrival[0] > self.missionWindow[1]:
                if verbose:
                    print("Mission window exceeded")
                break

            # Also break if the time of flight is too short (avoids singular
            # lambert solutions)
            tof = (
                time_at_arrival[i]
                - time_at_arrival[i - 1]
                - time_spent_preparing[i - 1]
            )
            if tof < 0.1:
                if verbose:
                    print("Time of flight too short or reached of chain.")
                break

            # Compute the ephemeris of the asteroid we are departing
            r1, v1 = self.asteroids[previous_ast_id].eph(
                T_START.mjd2000 +\
                time_at_arrival[i - 1] +\
                time_spent_preparing[i - 1]
            )

            # Compute the ephemeris of the next target asteroid
            r2, v2 = self.asteroids[current_ast_id].eph(
                T_START.mjd2000 + time_at_arrival[i]
            )

            # Solve the lambert problem for this flight
            l = pk.lambert_problem(
                r1=r1,
                r2=r2,
                tof=tof * pk.DAY2SEC,
                mu=self.MU,
                cw=False,
                max_revs=0)

            # Compute the delta-v necessary to go there and match its velocity
            DV1 = [a - b for a, b in zip(v1, l.get_v1()[0])]
            DV2 = [a - b for a, b in zip(v2, l.get_v2()[0])]
            DV = np.linalg.norm(DV1) + np.linalg.norm(DV2)

            # Compute propellant used for this transfer and update ship
            # propellant level
            propellant = propellant - DV / DV_PER_PROPELLANT

            # Break if we ran out of propellant during this transfer
            if propellant < 0:
                if verbose:
                    print("Out of propellant")
                break

            ###################### Step 2 #######################
            # If we are here, this asteroid-to-asteroid         #
            # jump is possible and we accumulate the mining     #
            # resource to the objective function.               #
            #####################################################

            # Get material of the asteroid we are visiting
            mat_idx = self.asteroidMaterialsTypes[current_ast_id]

            # Prepare as much material as is there or we have time to
            material_prepared[mat_idx] += np.minimum(
                self.asteroidMasses[current_ast_id],
                time_spent_preparing[i] / TIME_TO_MINE_FULLY,
            )
            material_prepared_at_t.append(material_prepared.copy())

            # If this is a propellant asteroid, we add it to the propellant
            if mat_idx == 3:
                propellant_found = np.minimum(
                    self.asteroidMasses[current_ast_id],
                    time_spent_preparing[i] / TIME_TO_MINE_FULLY,
                )
                propellant = np.minimum(1.0, propellant + propellant_found)

            if verbose:
                tank = f"{material_prepared[0]:.2f}|{material_prepared[1]:.2f}|\
                    {material_prepared[2]:.2f}"
                score = np.min(material_prepared[:3])
                print(
                    f"{current_ast_id}\t{time_at_arrival[i]:<4.2f}\
                    \t{propellant:<14.2f}\t{DV:<8.2f}\t{mat_idx}\t {tank}\
                    \t\t{score:.2f}"
                )

            visited = visited + 1
            propellant_at_time.append(propellant)

        return (
            material_prepared,
            material_prepared_at_t,
            ast_ids,
            time_at_arrival,
            time_spent_preparing,
            visited,
            propellant_at_time,
        )

    def fitness(self, x, verbose=False):
        """Evaluate the fitness of the decision variables.

        Args:
            x (numpy.array): Chromosome for the decision variables.
            verbose (bool): If True, print some info.

        Returns:
            float: Fitness of the chromosome.
        """

        # Evaluate the journey (checks feasibility and gets the material
        # prepared etc.)
        (
            material_prepared,
            _,
            ast_ids,
            time_at_arrival,
            time_spent_preparing,
            visited,
            _,
        ) = self._evaluate_journey(x, verbose)

        # The objective function in the end is the minimum
        # prepared mass of the three non-propellant material types.
        obj = np.min(material_prepared[:3])

        # Now the constraints
        # The visited asteroid ids must all be different (equality)
        ec = len(set(ast_ids[:visited])) - len(ast_ids[:visited])
        # The visiting epoch must be after the previous visiting epoch plus
        # the mining time (inequalities)
        ic = [0] * ((len(x) // 3) - 1)
        for i in range(1, visited):
            ic[i] = (
                time_at_arrival[i - 1]
                + time_spent_preparing[i - 1]
                - time_at_arrival[i]
            )
        return [-obj] + [ec] + ic

    def pretty(self, chromosome):
        """Pretty print the chromosome.

        Args:
            x (numpy.array): Chromosome for the decision variables.

        Returns:
            str: Pretty print of the chromosome.
        """
        self.fitness(chromosome, True)

    def example(self):
        """Returns an example solution."""
        # (disable automatic formatting for this)
        # fmt: off
        t_arr = [0,
                11.0,
                45.98091676982585,
                98.86574387748259,
                144.3421379448264,
                178.78720680368133,
                198.49061810149578,
                236.39180345018394,
                268.4772894184571]
        t_m = [0,
               18.980916769828053,
               22.88482710766111,
               29.47639406736512,
               17.445068858837555,
               18.703411297804774,
               19.901185348707877,
               24.085485968277332,
               17.543366859589646]
        a = [0,
             1446,
             5131,
             4449,
             8091,
             1516,
             151,
             4905,
             8490]
        # fmt: on
        return convert_to_chromosome(t_arr + t_m + a)

    def plot(self, solution):
        """Plots the journey of the solution and prepared materials
           and propellant.

        Args:
            solution (numpy.array): Chromosome for the decision variables.
        """
        # Evaluate the journey
        (
            _,
            material_prepared_at_t,
            _,
            time_at_arrival,
            time_spent_preparing,
            visited,
            propellant_at_time,
        ) = self._evaluate_journey(solution, verbose=True)

        # Throw away the first asteroid since it does not count
        time_at_arrival = time_at_arrival[1:]
        time_spent_preparing = time_spent_preparing[1:]

        # Plot the journey
        plt.figure(figsize=(10, 4), dpi=150)
        obj_function_scores = []
        for idx in range(visited):
            obj_function_scores.append(np.min(material_prepared_at_t[idx]))
        plt.plot(time_at_arrival[:visited], obj_function_scores, "--")

        lw = 8 * 8 // visited
        for idx in range(visited):
            # fmt:off
            plt.vlines(time_at_arrival[idx],0,material_prepared_at_t[idx][0],
                color="red",lw=lw)
            plt.vlines(time_at_arrival[idx] +\
                 3,0,material_prepared_at_t[idx][1],
                color="green",lw=lw)
            plt.vlines(time_at_arrival[idx] +\
                 6,0,material_prepared_at_t[idx][2],
                color="blue",lw=lw)
            plt.vlines(time_at_arrival[idx] - 3,0,propellant_at_time[idx],
                color="black",lw=lw)
            plt.hlines(-0.1,time_at_arrival[idx],time_at_arrival[idx] +\
                 time_spent_preparing[idx],color="grey",lw=6)

        plt.legend(
            [
                "Score",
                "Material 1",
                "Material 2",
                "Material 3",
                "Propellant",
                "Time preparing",
            ],
            bbox_to_anchor=(1.025, 0.5),
            loc="center left",
        )
        plt.xlabel("Time (days)")
        plt.tight_layout()
        plt.show()


udp = BeltMiningUdp([0, 1827.0])
result = udp.example()
udp.pretty(result)
