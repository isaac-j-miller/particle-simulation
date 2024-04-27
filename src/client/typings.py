from dataclasses import dataclass

from ..creator.typings import CreationConfig
from ..view.typings import ViewConfig

@dataclass
class ClientConfig:
    simulation_time_per_frame: float
    frames_per_second: int
    simulation_duration_seconds: int


@dataclass
class GlobalConfig:
    creation_config: CreationConfig
    client_config: ClientConfig
    view_config: ViewConfig
    gravitational_constant: float=10
    enable_collisions: bool=True
    
