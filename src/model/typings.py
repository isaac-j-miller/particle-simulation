from dataclasses import dataclass

@dataclass
class Particle:
    position: tuple[float, float]
    velocity: tuple[float, float]
    radius: float = 1
    mass: float = 1
    