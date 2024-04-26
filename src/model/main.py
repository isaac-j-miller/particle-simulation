from dataclasses import dataclass
from typing import Union
import math
from numpy import linalg as la, array

from .typings import Particle

FLOAT_ZERO = 1e-12

FLOAT_TOLERANCE = 1e-6

def float_equals(a: float, b: float) -> bool:
    return abs(a-b) <= FLOAT_TOLERANCE

@dataclass
class PreliminaryVector:
    current_position: list[float, float]
    next_position: list[float, float]
    difference_vector: list[float, float]

# ((vx, x0), x1), ((vy, y0), y1)
def to_equation(vect: PreliminaryVector) -> tuple[tuple[tuple[float, float], float]]:
   x_parametric = ((vect.difference_vector[0], -1), -1*vect.current_position[0])
   y_parametric = ((vect.difference_vector[1], -1), -1*vect.current_position[1])
   return x_parametric, y_parametric

def solve_intersection(eq1: tuple[tuple[float, float], float], eq2: tuple[tuple[float, float], float]) -> Union[None, float]:
    left_1, right_1 = eq1
    left_2, right_2 = eq2
    left =  array([left_1, left_2])
    right = array([right_1, right_2])
    try:
        return la.solve(left, right)
    except la.LinAlgError:
        return None
    

    # v1', v2'
    # [additive_solution, subtractive_solution]
def solve_momentum_equation(m1: float, m2: float, v1: float, v2: float) -> tuple[tuple[float, float], tuple[float, float]]:
    p_0 = m1*v1+m2*v2
    # p_0 = m1v1' + m2v2'
    k_0 = 0.5*m1*v1**2 + 0.5*m2*v2**2
    # k_0 = 0.5m1v1'**2 + 0.5m2v2'**2
    v2_prime_denominator = m2*(m2+m1)
    v2_prime_numerator_1st = p_0*m2
    v2_prime_numerator_under_sqrt = m1*m2*(2*k_0*m2+2*k_0*m1 - p_0**2)
    add_soln_v2_prime = (v2_prime_numerator_1st + math.sqrt(v2_prime_numerator_under_sqrt))/v2_prime_denominator
    sub_soln_v2_prime = (v2_prime_numerator_1st + math.sqrt(v2_prime_numerator_under_sqrt))/v2_prime_denominator
    add_soln_v1_prime = (p_0 - m2*add_soln_v2_prime)/m1
    sub_soln_v1_prime = (p_0 - m2*sub_soln_v2_prime)/m2
    return (add_soln_v1_prime, add_soln_v2_prime), (sub_soln_v1_prime, sub_soln_v2_prime)
    
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
        # 
        preliminary_vects: list[PreliminaryVector] = []

        for f_x, f_y, p in zip(forces_x, forces_y, self.particles):
            accel_x = f_x/p.mass
            accel_y = f_y/p.mass
            delta_v_x = accel_x * time_delta
            delta_v_y = accel_y * time_delta
            v_x = p.velocity[0] + delta_v_x
            v_y = p.velocity[1] + delta_v_y
            p_x = p.position[0] + v_x * time_delta
            p_y = p.position[1] + v_y * time_delta
            preliminary = PreliminaryVector(p.position, [p_x, p_y], [v_x, v_y])
            preliminary_vects.append(preliminary)
        
        done_collisions = set()
        for i, p1 in enumerate(preliminary_vects):
            particle_1 = self.particles[i]
            # TODO: allow for > 2 collisions
            soonest_collision_coords = None
            soonest_collision_time = None
            soonest_collision_idx = None
            for j, p2 in enumerate(preliminary_vects):
                if i == j:
                    continue
                hash_1 = f"{i}_{j}"
                hash_2 = f"{j}_{i}"
                if hash_1 in done_collisions or hash_2 in done:
                    continue
                p1_x, p1_y = to_equation(p1)
                p2_x, p2_y = to_equation(p2)
                intersection_x = solve_intersection(p1_x, p2_x)
                if intersection_x is not None and intersection_x[0] >=-1*FLOAT_TOLERANCE and intersection_x[0] < time_delta + FLOAT_TOLERANCE:
                    intersection_y = solve_intersection(p1_y, p2_y)
                    if (
                        intersection_y is not None and intersection_y[0] >=0 and float_equals(intersection_y[0],intersection_x[0])
                        ) and (
                            soonest_collision_time is None or soonest_collision_time > intersection_x[0]
                            ):                        
                        soonest_collision_time=intersection_x[0]
                        soonest_collision_coords=[intersection_x[1], intersection_y[1]]
                        soonest_collision_idx = j

            if soonest_collision_time is None or soonest_collision_coords is None or soonest_collision_idx is None:
                done_collisions.add(hash_1)
                particle_1.position = p1.next_position
                particle_1.velocity = p1.difference_vector
                continue
            particle_2 = self.particles[soonest_collision_idx]
            v_x_add, v_x_sub = solve_momentum_equation(particle_1.mass, particle_2.mass, p1.difference_vector[0], p2.difference_vector[0])
            v_y_add, v_y_sub = solve_momentum_equation(particle_1.mass, particle_2.mass, p1.difference_vector[1], p2.difference_vector[1])
            # TODO: figure out why we have 2 solutions, and which one to use. 
            # for now, use the additive solution
            # print(f"X add: ({v_x_add[0]}, {v_x_add[1]}); sub: ({v_x_sub[0]}, {v_x_sub[1]});")
            # print(f"Y add: ({v_y_add[0]}, {v_y_add[1]}); sub: ({v_y_sub[0]}, {v_y_sub[1]});")
            remaining_time = time_delta - soonest_collision_time
            p1_x = v_x_sub[0] * remaining_time + soonest_collision_coords[0]
            p1_y = v_y_sub[0] * remaining_time + soonest_collision_coords[1]
            particle_1.position = [p1_x, p1_y]
            particle_1.velocity = [v_x_sub[0], v_y_sub[0]]

            p2_x = v_x_sub[1] * remaining_time + soonest_collision_coords[0]
            p2_y = v_y_sub[1] * remaining_time + soonest_collision_coords[1]
            particle_2.position = [p2_x, p2_y]
            particle_2.velocity = [v_x_sub[1], v_y_sub[1]]
            done_collisions.add(hash_1)
        
        
