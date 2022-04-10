# SpockMining

![6cbf5b2f6921469d8d8a28f29345a248 large](https://user-images.githubusercontent.com/55914877/162574478-3efc0658-084e-49a7-9656-7d0fc39585f5.jpeg)


## Summary
The _Advanced Concepts Team_ (ACT) from _European Space Agency_ (ESA) partners with _Genetic and Evolutionary Computation Conference_ (GECCO) to propose a set of competitions.

These competitions are based into a future where the humanity, of several sophisticated robotic probes, has established a human settlement in another solar system: TRAPIST-1.

## Mining the belt
During the development of human settlement, a asteroid belt has been discoverd beyond the orbit of the outermost planet. These asteroids contain a rich source of materials considered essentials for setlling down including the material used for the probe propulsion itself.

### The challenge

_Mining the belt_ is single-objective optimization challenge with the goal of using your ship to prepare the material of as many as possible of 10000 asteroids into the belt outside the orbits of planets in the TRAPPIST-1 system for mining. There are three differents materials __Gold__,  __Platinum__ and  __Nickel__, and also the material required for refuelling the ship.

This mission has a deadling of __1827__ days to visit as many asteroid as you can.

Each asteroid provided a mass between __0.1__ or __1.0__ of the three materials of propellant.

The ship can store propellant for a delta-V of 10.000 m/s.

For each day at the asteroid an absolut __1/30__ mass is prepared/ propellant collected.

#### Constraints

This challenge has a set of constraints:

1. Each asteroid can only be visited once
2. Arrival times for visited asteroids need to be chronological
3. You can stay at one asteroid between 1 and 60 days
4. You have 1827 days to prepare
5. If you run out of propellant during a transfer to an asteroid, all following transfers are neglected.

#### Objective
The objective of this challenge is to maximize the minimum of the accumulated materials from all prepared asteroids (__gold__, __platinum__ and __nickel__)
