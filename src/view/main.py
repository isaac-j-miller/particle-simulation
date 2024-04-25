import pygame as pg

from ..model.typings import Particle
from ..view.typings import ViewConfig

class View:
    config: ViewConfig
    screen: pg.Surface
    heaviest_particle_idx: int
    def __init__(self, config: ViewConfig, heaviest_particle_idx: int):
        self.config = config
        self.heaviest_particle_idx = heaviest_particle_idx
        pg.init()
        pg.display.init()
        self.screen =pg.display.set_mode(self.config.dimensions)
        self.screen.fill([0,0,0])
        

    def display(self, particles: list[Particle]):
        self.screen.fill([0,0,0])
        transform = [0, 0]
        scaled_transform = [int(self.config.dimensions[0]/2), int(self.config.dimensions[1]/2)]
        if self.config.center_heaviest:
            heaviest = particles[self.heaviest_particle_idx].position
            transform = [-1*heaviest[0], -1*heaviest[1]]
        
        for p in particles:
            scaled = [int(((p.position[0]+transform[0])/ self.config.scale)) + scaled_transform[0], int((p.position[1]+transform[1])/self.config.scale)+scaled_transform[1]]
            self.screen.set_at(scaled, [255, 255, 255])
        pg.display.update()
