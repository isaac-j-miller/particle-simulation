from dataclasses import dataclass
from typing import Union


@dataclass
class ViewConfig:
    dimensions: tuple[int, int]
    scale: float
    center_heaviest: bool
    exile_distance: Union[None, float]