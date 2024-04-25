from dataclasses import dataclass


@dataclass
class ViewConfig:
    dimensions: tuple[int, int]
    scale: float
    center_heaviest: bool