import argparse
import jsons
import json

from src.client.main import Client
from src.client.typings import GlobalConfig
from src.creator.main import ParticleCreator
from src.model.main import SimulationModel
from src.view.main import View


def main(filename: str):
    with open(filename, "r") as f:
        content = json.load(f)
        config = jsons.load(content, GlobalConfig)
    creator = ParticleCreator(config.creation_config)
    particles = creator.create_initial_particles()
    max_mass = 0
    heaviest_particle_idx = 0
    for i, p in enumerate(particles):
        if p.mass >= max_mass:
            heaviest_particle_idx = i
            max_mass = p.mass
    model = SimulationModel(particles, config.gravitational_constant, config.repulsive_constant)
    view = View(config.view_config, heaviest_particle_idx)
    client = Client(view, model, config.client_config)
    view.display(particles)
    # input("Press enter to start")
    client.start()
    input("Press enter to exit")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file')
    args = parser.parse_args()
    main(args.config_file)
    exit(0)