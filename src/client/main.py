import pygame as pg
from .typings import ClientConfig
from ..model.main import SimulationModel
from ..view.main import View


class Client:
    model: SimulationModel
    view: View
    config: ClientConfig
    def __init__(self, view: View, model: SimulationModel, config: ClientConfig):
        self.view = view
        self.model = model
        self.config = config

    def start(self):
        clock = pg.time.Clock()
        total_frames = self.config.simulation_duration_seconds * self.config.frames_per_second
        for _ in range(total_frames):
            self.model.calculate_new_positions(self.config.simulation_time_per_frame)
            to_exile = self.view.display(self.model.particles, self.model.exiled_particles)
            for idx in to_exile:
                self.model.exile_particle(idx)
            clock.tick(self.config.frames_per_second)




