from typing import Union
import pygame as pg

from numpy import array
from src.model.main import magnitude

from ..model.typings import Particle
from ..view.typings import ViewConfig

class View:
    config: ViewConfig
    screen: pg.Surface
    heaviest_particle_idx: int
    exile_distance: Union[float, None]
    def __init__(self, config: ViewConfig, heaviest_particle_idx: int):
        self.config = config
        self.heaviest_particle_idx = heaviest_particle_idx
        pg.init()
        pg.display.init()
        self.screen =pg.display.set_mode(self.config.dimensions)
        self.screen.fill([0,0,0])
        self.exile_distance = magnitude(self.config.dimensions) * self.config.scale / 2 + self.config.exile_distance if self.config.exile_distance is not None else None
        

    def display(self, particles: list[Particle], exiled: set[int]) -> set[int]:
        self.screen.fill([0,0,0])
        transform = [0, 0]
        scaled_transform = [int(self.config.dimensions[0]/2), int(self.config.dimensions[1]/2)]
        center_to_use = [int(self.config.dimensions[0]/2), int(self.config.dimensions[1]/2)]
        if self.config.center_heaviest:
            heaviest = particles[self.heaviest_particle_idx].position
            center_to_use = heaviest
            transform = [-1*heaviest[0], -1*heaviest[1]]
        to_exile = set()
        for i, p in enumerate(particles):
            if i in exiled:
                continue
            scaled = [int(((p.position[0]+transform[0])/ self.config.scale)) + scaled_transform[0], int((p.position[1]+transform[1])/self.config.scale)+scaled_transform[1]]
            color = [255, 255, 255] if i % 2 == 0 else [255,255,0] 
            if self.exile_distance is not None:
                distance_from_center = magnitude(array(p.position)-array(center_to_use))
                if distance_from_center > self.exile_distance:
                    print(f"exiling {i} because distance from center of view > {int(self.exile_distance)}")
                    to_exile.add(i)
                    continue
            try:
                scaled_rad = int(p.radius/self.config.scale)
                if(scaled_rad <=1):
                    self.screen.set_at(scaled, color)
                else:
                    pg.draw.circle(self.screen, color, scaled,scaled_rad)
            except OverflowError:
                continue
        pg.display.update()
        return to_exile
