import argparse
import sys
from typing import Union
import jsons
import json

from src.client.main import Client
from src.client.typings import GlobalConfig
from src.creator.main import ParticleCreator
from src.model.main import SimulationModel
from src.view.main import View


def main(filename: str, output_file: Union[str, None]):
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
    model = SimulationModel(particles, config.gravitational_constant, config.enable_collisions)
    view = View(config.view_config, heaviest_particle_idx)
    client = Client(view, model, config.client_config)
    view.display(particles, set())
    client.start()
    print("Simulation complete")
    if output_file is None:
        save_confirmed = None
        while save_confirmed is None:
            should_save = input("Would you like to save this as an image? (y/n): ")
            if should_save == "y":
                save_confirmed = True
            elif should_save == "n":
                save_confirmed = False
            else:
                print("Invalid input")
        if save_confirmed:
            while True:
                f = input("Please input the filepath to save the image to: ")
                try:
                    view.save_image(f)
                    return
                except Exception:
                    the_type, the_value, _the_traceback = sys.exc_info()
                    print(f"Error saving file: {the_type} {the_value}")
    else:
        view.save_image(output_file)

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file')
    parser.add_argument('-o', '--image-output', default=None)
    args = parser.parse_args()
    main(args.config_file, args.image_output)
    exit(0)