from typing import Literal, List
from dataclasses import dataclass

from ..model.typings import Particle


@dataclass
class DistributionParams:
    mode: str

@dataclass
class VectorParams:
    random: bool
    range: tuple[float, float]
    base_value: tuple[float, float]
    
@dataclass
class MassParams:
    random: bool
    range: tuple[float, float]
    base_value: float
    
@dataclass
class RandomParticleDistributionParams(DistributionParams):
    mode: Literal["random"]
    x_range: tuple[float, float]
    y_range: tuple[float, float]

@dataclass
class LinearParticleDistributionParams(DistributionParams):
    mode: Literal["linear"]
    angle: float
    starting_point: tuple[float, float]
    spacing: float

@dataclass
class CircularParticleDistributionParams(DistributionParams):
    mode: Literal["circular"]
    fill: bool
    radius: float
    center_point: tuple[float, float]

@dataclass
class RadiusParams:
    random: bool
    range: tuple[float, float]
    base_value: float

@dataclass
class GenerationParams:
    num_particles: int
    distribution_params: DistributionParams
    vector_params: VectorParams
    mass_params: MassParams
    radius_params: RadiusParams


@dataclass
class CreationConfig:
    particle_generation: List[GenerationParams]
    known_particles: List[Particle]