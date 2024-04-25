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
    gravitational_constant: float
    repulsive_constant: float
    creation_config: CreationConfig
    client_config: ClientConfig
    view_config: ViewConfig
    
