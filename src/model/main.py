from dataclasses import dataclass
from typing import Union
import math
import numpy as np
from numpy import linalg as la, array

from .typings import Particle

@dataclass
class PreliminaryVector:
    current_position: list[float, float]
    next_position: list[float, float]
    difference_vector: list[float, float]

FLOAT_ZERO = 1e-12

FLOAT_TOLERANCE = 1e-6

def float_equals(a: float, b: float) -> bool:
    return abs(a-b) <= FLOAT_TOLERANCE

def eval_equation(v: tuple[float, float], c: tuple[float, float], t: float)-> tuple[float, float]:
    x = v[0]*t + c[0]
    y = v[1]*t + c[1]
    return (x, y)


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
    

def magnitude(vector: tuple[float, float]) -> float:
    return math.sqrt(vector[0]**2+vector[1]**2)

def quadratic_equation(a: float, b: float, c: float) -> Union[list[float,float], None]:
    if a == 0:
        return None
    discriminant = b**2-4*a*c
    if(discriminant < 0):
        return None
    if float_equals(discriminant, 0):
        return [-0.5*b/a]
    sqrt = math.sqrt(discriminant)
    add = (-1*b + sqrt)/(2*a)
    sub = (-1*b - sqrt)/(2*a)
    return [sub, add]

def get_collision_time(
        r1: float, r2: float, 
        vect1: PreliminaryVector,
        vect2: PreliminaryVector,
        time_delta: float) -> Union[None, float]:
    v1 = vect1.difference_vector
    v2 = vect2.difference_vector
    i1 = vect1.current_position
    i2 = vect2.current_position
    # create a quadratic equation to solve for t(d) where d = r1+r2
    # ax^2 + bx + c = 0
    a = v2[0]**2 - 2*v2[0]*v1[0] + v1[0]**2 + v2[1]**2 - 2*v2[1]*v1[1] + v1[1]**2
    b = 2*(v2[0]*i2[0] - i2[0]*v1[0] - v2[0]*i1[0] + v1[0]*i1[0] + v2[1]*i2[1] - i2[1]*v1[1] - v2[1]*i1[1] + v1[1]*i1[1])
    c = i2[0]**2 - 2 * (i2[0]*i1[0]+i2[1]*i1[1]) + i1[0]**2 + i2[1]**2 + i1[1]**2 - (r1+r2)**2
    intersection_times = quadratic_equation(a, b, c)
    if intersection_times is None:
        return None
    acceptable = [t for t in intersection_times if t >= -1*FLOAT_ZERO and t <= time_delta+FLOAT_ZERO]
    if len(acceptable) == 0:
        # they have already collided or won't collide within the time frame
        return None
    return acceptable[0]

def solve_2d_momentum_equation_2(
        m1: float, m2: float, 
        vect1: PreliminaryVector,
        vect2: PreliminaryVector,
        time_delta: float,
        intersection_time: float
        ) -> Union[None, tuple[PreliminaryVector, PreliminaryVector]]:
    v1 = vect1.difference_vector
    v2 = vect2.difference_vector
    i1 = vect1.current_position
    i2 = vect2.current_position
    c1 = eval_equation(v1, i1, intersection_time)
    c2 = eval_equation(v2, i2, intersection_time)
    v1_vect = array(v1)
    v2_vect = array(v2)
    c1_vect = array(c1)
    c2_vect = array(c2)
    fv1 = v1_vect - (2*m2/(m1+m2))*(np.dot(v1_vect-v2_vect, c1_vect-c2_vect)/magnitude(c1_vect-c2_vect)**2)*(c1_vect-c2_vect)
    fv2 = v2_vect - (2*m1/(m1+m2))*(np.dot(v2_vect-v1_vect, c2_vect-c1_vect)/magnitude(c2_vect-c1_vect)**2)*(c2_vect-c1_vect)
    dt_after_intersection=time_delta-intersection_time
    fp1 = eval_equation(fv1,c1, dt_after_intersection)
    fp2 = eval_equation(fv2,c2, dt_after_intersection)
    f1 = PreliminaryVector([c1[0], c1[1]], fp1, [fv1[0], fv1[1]])
    f2 = PreliminaryVector([c2[0], c2[1]], fp2, [fv2[0], fv2[1]])
    return f1, f2

class SimulationModel:
    particles: list[Particle]
    gravitational_constant: float
    def __init__(self, particles: list[Particle], gravitational_constant: float):
        self.particles = particles
        self.gravitational_constant = gravitational_constant

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
                theta = math.atan2(delta_y, delta_x)
                grav_force_x = gravitational_force * math.cos(theta)
                grav_force_y = gravitational_force * math.sin(theta)
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
        particles_processed = set()
        for i, p1 in enumerate(preliminary_vects):
            particle_1 = self.particles[i]
            # TODO: allow for > 2 collisions
            soonest_collision_time = None
            soonest_collision_idx = None
            for j, p2 in enumerate(preliminary_vects):
                if i == j:
                    continue
                hash_1 = f"{i}_{j}"
                hash_2 = f"{j}_{i}"
                if hash_1 in done_collisions or hash_2 in done:
                    continue
                intersection_time= get_collision_time(particle_1.radius, self.particles[j].radius, p1, p2, time_delta)
                if intersection_time is None:
                    done_collisions.add(hash_1)
                    continue
                soonest_collision_time = intersection_time
                soonest_collision_idx = j
                # print(f"imminent collision between {i},{j} at t={soonest_collision_time}")

            done_collisions.add(hash_1)
            if soonest_collision_time is None or soonest_collision_idx is None:
                continue
            particle_2 = self.particles[soonest_collision_idx]
            result = solve_2d_momentum_equation_2(particle_1.mass, particle_2.mass, p1, p2, time_delta, soonest_collision_time)
            if result is None:
                # print(f"invalid solution for collision of {i}, {soonest_collision_idx}")
                continue
            # print(f"valid solution for collision of {i}, {soonest_collision_idx}")
            p1_new, p2_new = result
            particle_1.position = p1_new.next_position
            particle_1.velocity = p1_new.difference_vector
            particle_2.position = p2_new.next_position
            particle_2.velocity = p2_new.difference_vector
            particles_processed.add(i)
            particles_processed.add(j)

        for i, prelim in enumerate(preliminary_vects):
            if i in particles_processed:
                continue
            p = self.particles[i]
            p.position = prelim.next_position
            p.velocity = prelim.difference_vector
        
        
