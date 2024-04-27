from dataclasses import dataclass
from typing import Union


@dataclass
class ViewConfig:
    dimensions: tuple[int, int]=800,800
    scale: float=1
    center_heaviest: bool=False
    exile_distance: Union[None, float]=5000
    color_velocity: bool = True
    trace: bool = False
    color_scheme: str = "viridis"
    use_relative_color_scale: bool = True
    color_scale_velocity_bounds: Union[None, tuple[float, float]] = None
    background_color_rgb: tuple[int, int, int] = 0,0,0
    default_particle_color_rgb: tuple[int, int, int] = 255,255,255