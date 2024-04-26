from dataclasses import dataclass

@dataclass
class Particle:
    mass: float
    position: list[float, float]
    velocity: list[float, float]
    radius: float = 1
    