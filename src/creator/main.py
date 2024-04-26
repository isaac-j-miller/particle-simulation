import random
import math
from ..model.typings import Particle
from .typings import CreationConfig, RandomParticleDistributionParams, CircularParticleDistributionParams, LinearParticleDistributionParams, GenerationParams

phi = (1.0 + math.sqrt(5)) / 2.0
phi_sq = phi**2.0


class ParticleGenerator:
    config: GenerationParams
    def __init__(self, config: GenerationParams):
        self.config = config

    @staticmethod
    def _random_float(min: float, max: float) -> float:
        base = random.random()
        return (base - min)/(max - min)
    
    def _get_particle_mass(self) -> float:
        if self.config.mass_params.random:
            return ParticleGenerator._random_float(self.config.mass_params.range[0], self.config.mass_params.range[1]) + float(self.config.mass_params.base_value)
        return float(self.config.mass_params.base_value)

    def _get_particle_radius(self) -> float:
        if self.config.radius_params.random:
            return ParticleGenerator._random_float(self.config.radius_params.range[0], self.config.radius_params.range[1]) + float(self.config.radius_params.base_value)
        return float(self.config.radius_params.base_value)

    def _get_particle_vector(self) -> list[float, float]:
        if self.config.vector_params.random:
            return [
                float(self.config.vector_params.base_value[0]) + ParticleGenerator._random_float(self.config.vector_params.range[0], self.config.vector_params.range[1]),
                float(self.config.vector_params.base_value[1]) + ParticleGenerator._random_float(self.config.vector_params.range[0], self.config.vector_params.range[1])]
        return [float(self.config.vector_params.base_value[0]), float(self.config.vector_params.base_value[1])]

    def _create_random_particles(self, distribution_config: RandomParticleDistributionParams) -> list[Particle]:
        particles: list[Particle] = []
        for _ in range(self.config.num_particles):
            mass = self._get_particle_mass()
            vector = self._get_particle_vector()
            radius = self._get_particle_radius()
            position = [ParticleGenerator._random_float(distribution_config.x_range[0], distribution_config.x_range[1]),ParticleGenerator._random_float(distribution_config.y_range[0], distribution_config.y_range[1])]
            particle = Particle(mass, position, vector, radius)
            particles.append(particle)
        return particles

    def _create_circular_particles(self, distribution_config: CircularParticleDistributionParams) -> list[Particle]:
        positions: list[list[float, float]] = []
        if distribution_config.fill:
            spacing_constant = distribution_config.radius / math.sqrt(self.config.num_particles - 1)
            for i in range(self.config.num_particles):
                theta = 2 * math.pi * (i+1) / phi_sq
                r = spacing_constant * math.sqrt(i+1)
                x = r * math.cos(theta) + distribution_config.center_point[0]
                y = r * math.sin(theta) + distribution_config.center_point[1]
                positions.append([x, y])
        else:
            spacing = 2 * math.pi / self.config.num_particles
            for i in range(self.config.num_particles):
                theta = spacing * i
                r = distribution_config.radius
                x = r * math.cos(theta) + distribution_config.center_point[0]
                y = r * math.sin(theta) + distribution_config.center_point[1]
                positions.append([x, y])
        particles: list[Particle] = []
        for position in positions:
                mass = self._get_particle_mass()
                vector = self._get_particle_vector()    
                radius = self._get_particle_radius()
                particle = Particle(mass, position, vector, radius)
                particles.append(particle)

        return particles

    def _create_linear_particles(self, distribution_config: LinearParticleDistributionParams) -> list[Particle]:
        particles: list[Particle] = []
        for i in range(self.config.num_particles):
            x = distribution_config.starting_point[0] + math.cos(distribution_config.angle) * distribution_config.spacing * i
            y = distribution_config.starting_point[1] + math.sin(distribution_config.angle) * distribution_config.spacing * i
            position = [x, y]
            mass = self._get_particle_mass()
            vector = self._get_particle_vector()    
            radius = self._get_particle_radius()
            particle = Particle(mass, position, vector, radius)
            particles.append(particle)


    def create_distribution(self) -> list[Particle]:
        if self.config.num_particles == 0:
            return []
        if self.config.distribution_params.mode == "random":
            return self._create_random_particles(self.config.distribution_params)
        if self.config.distribution_params.mode == "circular":
            return self._create_circular_particles(self.config.distribution_params)
        if self.config.distribution_params.mode == "linear":
            return self._create_linear_particles(self.config.distribution_params)
        raise ValueError(f"distribution mode {self.config.distribution_params.mode} is invalid")

class ParticleCreator:
    config: CreationConfig
    def __init__(self, params: CreationConfig):
        self.config = params
        
    def create_initial_particles(self) -> list[Particle]:
        existing = self.config.known_particles
        arr = [*existing]
        for params in self.config.particle_generation:
            creator = ParticleGenerator(params)
            particles = creator.create_distribution()
            arr+=particles
        return arr