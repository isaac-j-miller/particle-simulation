import math
from .typings import Particle

FLOAT_ZERO = 0.0000000000001

class SimulationModel:
    particles: list[Particle]
    gravitational_constant: float
    repulsive_constant: float
    def __init__(self, particles: list[Particle], gravitational_constant: float, repulsive_constant: float):
        self.particles = particles
        self.gravitational_constant = gravitational_constant
        self.repulsive_constant = repulsive_constant


    def calculate_new_positions(self, time_delta: float):
        forces_x: list[float] = [0 for _ in self.particles]
        forces_y: list[float] = [0 for _ in self.particles]
        done = set()
        for i, p1 in enumerate(self.particles):
            for j, p2 in enumerate(self.particles):
                if i == j:
                    continue
                hash_1 = f"{i}_{j}"
                hash_2 = f"{j}_{i}"
                if hash_1 in done or hash_2 in done:
                    continue
                delta_y = p2.position[1] - p1.position[1]
                delta_x = p2.position[0] - p1.position[0]
                distance_squared = delta_x**2 + delta_y**2
                if distance_squared == 0:
                    continue
                gravitational_force = self.gravitational_constant * p1.mass * p2.mass / distance_squared
                repulsive_force = self.repulsive_constant / distance_squared**2
                net_force = gravitational_force - repulsive_force
                theta = math.atan2(delta_y, delta_x)
                grav_force_x = net_force * math.cos(theta)
                grav_force_y = net_force * math.sin(theta)
                forces_x[i] += grav_force_x
                forces_y[i] += grav_force_y
                forces_x[j] -= grav_force_x
                forces_y[j] -= grav_force_y
                done.add(hash_1)
        for f_x, f_y, p in zip(forces_x, forces_y, self.particles):
            accel_x = f_x/p.mass
            accel_y = f_y/p.mass
            delta_v_x = accel_x * time_delta
            delta_v_y = accel_y * time_delta
            p.velocity[0] += delta_v_x
            p.velocity[1] += delta_v_y
            # TODO: don't set this here once collisions are implemented
            p.position[0] += p.velocity[0] * time_delta
            p.position[1] += p.velocity[1] * time_delta
        # TODO: collisions
        
        
