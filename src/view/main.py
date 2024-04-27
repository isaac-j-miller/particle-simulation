from typing import Union
import pygame as pg
from os import path, makedirs
from matplotlib import colormaps
from matplotlib.colors import Colormap
from numpy import array
from src.model.main import magnitude

from ..model.typings import Particle
from ..view.typings import ViewConfig


class View:
    config: ViewConfig
    screen: pg.Surface
    heaviest_particle_idx: int
    exile_distance: Union[float, None]
    cmap: Colormap
    def __init__(self, config: ViewConfig, heaviest_particle_idx: int):
        self.config = config
        self.heaviest_particle_idx = heaviest_particle_idx
        if self.config.color_velocity and not self.config.use_relative_color_scale and not self.config.color_scale_velocity_bounds:
            raise Exception("Invalid Config: view_config.color_scale_velocity_bounds must be defined if view_config.use_relative_color_scale is False and view_config.color_velocity is True")
        pg.init()
        pg.display.init()
        self.cmap = colormaps[self.config.color_scheme]
        self.screen =pg.display.set_mode(self.config.dimensions)
        self.screen.fill(self.config.background_color_rgb)
        self.exile_distance = magnitude(self.config.dimensions) * self.config.scale / 2 + self.config.exile_distance if self.config.exile_distance is not None else None
        
    def save_image(self, fp: str):
        dirname = path.dirname(fp)
        if dirname:
            makedirs(dirname, exist_ok=True)
        pg.image.save(self.screen, fp)

    def display(self, particles: list[Particle], exiled: set[int]) -> set[int]:
        if not self.config.trace:
            self.screen.fill(self.config.background_color_rgb)
        transform = [0, 0]
        scaled_transform = [int(self.config.dimensions[0]/2), int(self.config.dimensions[1]/2)]
        center_to_use = [int(self.config.dimensions[0]/2), int(self.config.dimensions[1]/2)]
        if self.config.color_velocity:
            if self.config.use_relative_color_scale:
                velocities = [magnitude(p.velocity) for i, p in enumerate(particles) if i not in exiled]
                if len(velocities) == 0:
                    return set()
                min_vel = min(velocities)
                max_vel = max(velocities)
            else:
                min_vel = self.config.color_scale_velocity_bounds[0]
                max_vel = self.config.color_scale_velocity_bounds[1]
            if self.config.center_heaviest:
                point_of_reference_velocity = array(particles[self.heaviest_particle_idx].velocity)
            else:
                point_of_reference_velocity = array([0, 0])
        if self.config.center_heaviest:
            heaviest = particles[self.heaviest_particle_idx].position
            center_to_use = heaviest
            transform = [-1*heaviest[0], -1*heaviest[1]]

        to_exile = set()
        for i, p in enumerate(particles):
            if i in exiled:
                continue
            if self.exile_distance is not None:
                distance_from_center = magnitude(array(p.position)-array(center_to_use))
                if distance_from_center > self.exile_distance:
                    # print(f"yeet {i}")
                    to_exile.add(i)
                    continue
            scaled = [int(((p.position[0]+transform[0])/ self.config.scale)) + scaled_transform[0], int((p.position[1]+transform[1])/self.config.scale)+scaled_transform[1]]
            if self.config.color_velocity:
                speed=magnitude(array(p.velocity)-point_of_reference_velocity)
                vel_normalized = (speed - min_vel)/(max_vel-min_vel)
                # if not using relative scale, vel_normalized may fall outside the bounds, so we must clip it
                if vel_normalized <= 0:
                    vel_normalized = 0.000001
                if vel_normalized >= 1:
                    vel_normalized = 0.999999
                color = [int(x) for x in array(self.cmap(vel_normalized)[:-1])*255]
            else:
                color = self.config.default_particle_color_rgb
            try:
                scaled_rad = int(p.radius/self.config.scale)
                if(scaled_rad <=1):
                    self.screen.set_at(scaled, color)
                else:
                    pg.draw.circle(self.screen, color, scaled,scaled_rad)
            except OverflowError:
                # this happens when a particle is VERY far outside the bounds. this should never happen if exile_distance is enabled.
                continue
        pg.display.update()
        return to_exile
