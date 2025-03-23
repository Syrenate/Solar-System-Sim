import os
import json
import math
import time
import numpy as np

import matplotlib.path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation

from BodyClasses import CelestialBody, Planet

# Things left to do:
#   - Experiment 1; alternative integration methods


class Simulation:
    """Class to represent the matplotlib simulation space, handles UI elements and simulation rendering."""
    def __init__(self, reference_space, reference_object):
        """Constructor method for the simulation; defines the simulation and UI axes in matplotlib, according to a reference Space instance."""
        colour = 'white'
        plt.rcParams['text.color'] = colour
        plt.rcParams['axes.labelcolor'] = colour
        plt.rcParams['xtick.color'] = colour
        plt.rcParams['ytick.color'] = colour

        self.space = reference_space        # Keep reference to the instance of Space we're simulating.
        self.space_reference_object = reference_object
        self.sim_fig, self.sim_ax = plt.figure(0), self.generate_simulation_axes()   # Initialise the simulation axes.
        self.sim_annotations = []

        self.sim_fig.set_facecolor((0.12,0.12,0.12))
        self.sim_ax.patch.set_edgecolor('white')
        self.sim_ax.patch.set_linewidth(2)
        self.sim_ax.set_facecolor((0.02,0.02,0.02))
        plt.title("Many-Body Simulation of our Solar System", color='white')

        def update_simulation_limit(val):   # Method to rescale the simulation boundary to zoom in/out.
            scale = self.space.greatest_orbital_radius / val
            self.sim_ax.set_xlim(-scale, scale)
            self.sim_ax.set_ylim(-scale, scale)

        def update_time_scale(val):     # Method which updates the space's time scale.
            self.space.time_step = self.space.original_time_step * val

        self.space.current_energy_level = sum(self.space.get_energy_levels())

        self.sim_ax.axis("scaled")
        update_simulation_limit(5)
        self.sim_fig.subplots_adjust(left=0.3, bottom = 0.3)


        axzoom = self.sim_fig.add_axes([0.35, 0.15, 0.55, 0.03])
        self.zoom_slider = Slider(      # Interactive slider to zoom in/out of the simulation space.
            ax=axzoom,
            label='Zoom scale',
            valmin=0.5,
            valmax=10,
            valinit=5,
        )

        axtime = self.sim_fig.add_axes([0.35, 0.05, 0.55, 0.03], xscale='log')
        self.time_slider = Slider(      # Interactive time scale slider for the space instance.
            ax=axtime,
            label='Time scale',
            valmin=0.1,
            valmax=20,
            valinit=1,
        )

        self.zoom_slider.on_changed(update_simulation_limit)
        self.time_slider.on_changed(update_time_scale)

    def run(self):
        """Method called to initialise the simulation."""
        self.anim = FuncAnimation(self.sim_fig, lambda i: self.space.update(i, self), frames=self.space.iteration_limit, repeat=False, interval=1/3000, blit=True)
        plt.show()

    def update_patch_positions(self, reference):
        for i in range(len(self.space.bodies)):
            body = self.space.bodies[i]
            self.patches[i].center = (body.position - reference).get()
            if body.type == "planet": self.patches[i + self.space.body_count - 1].set_path(matplotlib.path.Path(body.position_trail))

    def generate_simulation_axes(self):
        """Method to create patches for each simulation element; all physical objects and their trails if applicable."""
        ax = plt.axes()
        self.patches = []
        for body in [self.space.star] + self.space.planets:     # Render each celestial body
            self.patches.append(mpatches.Circle(body.position.get(), 0.2 if body.name == "sun" else 0.1, color = body.colour, animated = True))
            ax.add_patch(self.patches[-1])

        for planet in self.space.planets:     # Render the trail of each planet
            planet.position_trail.append(planet.position.get())
            path = matplotlib.path.Path(planet.position_trail, codes = [np.uint8(1)], closed=True)
            self.patches.append(mpatches.PathPatch(path, color = planet.colour, animated = True, fill = False))
            ax.add_patch(self.patches[-1])

        return ax

    def generate_output(self):
        """Method which generates a console output of simulation and experiment data."""
        output = ""

        simulation_output =  f"     Elapsed time: {round(self.space.elapsed_time, 3)} Years (Current frame: {self.space.elapsed_frames}/{self.space.iteration_limit})\n"
        simulation_output += f"     Current framerate: {round(np.mean(self.space.framerate_history), 3)}fps, Frame speed: {round(self.space.time_step, 3)} Y/frame."
        output += "Simulation data:\n" + simulation_output + "\n\n"

        # Experiment 1
        planet_output = ""
        for planet in self.space.planets: planet_output += planet.orbit_output
        if len(planet_output) > 0: output += "Orbital Periods (Experiment 1):\n" + planet_output + "\n"

        # Experiment 2
        k_energy, p_energy = self.space.get_energy_levels()
        energy_output =  f"     Kinetic Energy: {round(k_energy, 3)}, Potential Energy: {round(p_energy, 3)}\n"
        energy_output += f"     Total Energy: {round(k_energy + p_energy, 3)}. "
        energy_output += f"Recent standard deviation [sample size: {len(self.space.energy_level_history)}]: {round(math.sqrt(np.var(self.space.energy_level_history)), 3)}"
        output += "Energy Levels (Experiment 2):\n" + energy_output + "\n\n"

        # Experiment 4
        allignmnet_output = f"Allignment occourances (Experiment 4): \n"
        allignmnet_output += f"     Years at which planets have alligned: {self.space.allignment_history}" 
        output += allignmnet_output

        return output

class Space:
    """Class in which the simulation is contained; all simulation and rendering occours here."""
    def __init__(self, solar_parameters, file):
        """Constructor method in which the solar parameters are parsed to initialise the many-body system."""
        self.iteration_limit = solar_parameters['num_iterations']   # Maximum iterations for the simulation. Terminates once reached.
        self.time_step = solar_parameters['timestep']               # Time step; the time interval to be simulated each frame.
        self.original_time_step = solar_parameters['timestep']      # Copy of the time_step parameter for use once the main value is changed.
        self.grav_constant = solar_parameters['grav_const']
        self.default_os = solar_parameters['default_os']
        self.integration_type = solar_parameters['integration_method']

        self.greatest_orbital_radius = 0    # Records greatest orbit to set axis limits in the simulaiton.
        self.json_file = file               # Reference to the json file so that data may be written in.

        self.elapsed_time = 0               
        self.elapsed_frames = 0                  
        self.framerate_history = []
        self.framerate_history_length = 30

        self.energy_level_history = []           # Record of previous energy levels for study.
        self.energy_history_sample_size = 1000   # Desired maximum length of energy_level_history.
        self.energy_sample = 100                 # When recording the energy level, only sample every energy_sample frames.
        self.trail_sample = 3                    # Same concept as above for sampling the position trail.

        self.allignment_history = []      # Record occourances of planetary allignment.
        self.allignment_sample = solar_parameters['experiment_4']['planet_sample_size']        # How many planets should be considered for allignment.
        self.allignment_tolerance = solar_parameters['experiment_4']['allignment_tolerance'] * math.pi / 180    # Greatest angle between a planet and the mean angle to be considered alligned.
 
        self.planets = []       # Array of planets in the system.
        self.body_count = 0     # Total number of celestial bodies.
        for body in solar_parameters['bodies']:     # This loop scans all bodies specified in the json file, and creates instances to represent them.
            if body['name'] == "sun": self.star = CelestialBody(self, body, "star")    # self.star will reference the sun in our solar system.
            else: self.planets.append(Planet(self, body, body['nasa_orbit_period']))
            self.body_count += 1 

            # Compare radius to the greatest radius recorded; replace if greater. This ensures all planets are on screen at 100% zoom.
            if body['orbital_radius'] > self.greatest_orbital_radius: self.greatest_orbital_radius = body['orbital_radius']

        for planet in self.planets:     # Once the greatest radius has been evaluated, set the trail length 
            planet.position_trail_length = math.ceil(3000 * (planet.position.x / self.greatest_orbital_radius))
            
        self.bodies = [self.star] + self.planets    # Array of all celestial bodies

    def update(self, iteration_count, simulation):
        start_update_time = time.time()

        self.elapsed_time += self.time_step
        self.elapsed_frames += 1
        """Method to simulate planetary motion in one time step."""
        self.simulation = simulation
        self.determine_allignment(self.allignment_sample)
        bodies_copy = self.bodies.copy()

        for i, body in enumerate(self.bodies):   # This loop updates all bodies in space, and evaluates if a planet has orbited or not.
            prev_pos = body.position
            if self.integration_type == 'beeman_integration': body.beeman_update(i == len(self.bodies) - 1)
            elif self.integration_type == 'euler_cromer': body.euler_cromer_update(bodies_copy)
            elif self.integration_type == 'direct_euler': body.direct_euler_update(bodies_copy)

            if body.type == "planet":
                # If not the first frame, the planet has not revolved before, the planet is in the +ve x-axis region, and it crosses the axis boundary; it has revolved. Record orbital period.
                if iteration_count != 0 and not body.has_revolved and body.position.x > 0 and (body.position.y > 0 and prev_pos.y <= 0):
                    body.orbital_period = self.elapsed_time    # Orbital period is the time since the simulation started.
                    body.has_revolved = True    
                    body.gen_orbit_output()

                body.angle_from_xaxis = body.position.angle_between_xaxis()     # Simplified from dot product formula with (1,0), the +ve x-axis.

            if body.type != "star":     # Every non-star body should have a trail (since the star should be centered).
                if iteration_count % self.trail_sample == 0: body.position_trail.append((body.position - self.star.position).get())
                if len(body.position_trail) > body.position_trail_length: body.position_trail = body.position_trail[1:]

        energy = sum(self.get_energy_levels())    # Evaluate current energy level of the system & add to history.
        if iteration_count % self.energy_sample: 
            self.energy_level_history.append(energy)
            if len(self.energy_level_history) > self.energy_history_sample_size: 
                self.energy_level_history = self.energy_level_history[1:]

        self.simulation.update_patch_positions(self.star.position)
        time_delay = time.time() - start_update_time
        if time_delay != 0: self.framerate_history.append(1/time_delay)
        if len(self.framerate_history) > self.framerate_history_length: self.framerate_history = self.framerate_history[1:]


        if iteration_count % 30 == 0:
            if self.default_os == "windows": os.system('cls')
            elif self.default_os == "linux": os.system('clear')
            print(self.simulation.generate_output())
        return simulation.patches
    
    def determine_allignment(self, sample_size):
        """Method which checks planetary allignment"""
        allignments = []
        for planet in self.planets[:sample_size]:    # First generate the angle from 0 to 2pi of each planet.
            angle = planet.position.angle_between_xaxis()
            if planet.position.y < 0: angle = 2 * math.pi - angle
            allignments.append(angle)

        mean_angle = np.mean(allignments)
        alligned = True     # Now loop through each plannet to check if they are close enough to the mean angle. Is so, they are 'alligned'.
        for i, planet in enumerate(self.planets[:sample_size]):
            if alligned:
                angle = allignments[i]
                angle_min, angle_max = min([angle, mean_angle]), max([angle, mean_angle])
                min_angle_delta = min([angle_max - angle_min, 2 * math.pi - (angle_max - angle_min)])
                if min_angle_delta > self.allignment_tolerance: alligned = False

        if len(self.allignment_history) == 0 or (alligned and np.abs(self.allignment_history[-1] - round(self.elapsed_time,2)) > 0.2): 
            self.allignment_history.append(float(round(self.elapsed_time, 2)))


    def get_energy_levels(self):
        """Method that calculates and returns the kinetic & potential energy levels of all bodies."""
        k_energy, p_energy = 0, 0
        for body in self.bodies:
            k_energy += body.mass * (body.velocity.magnitude() ** 2)
            p_energy -= 0 if body.type == "star" else body.mass / body.position.magnitude()
        return k_energy / 2, p_energy * (self.star.mass * self.grav_constant)
            

if __name__ == "__main__":
    with open("simulation_parameters.json") as param_file:
        space = Space(json.load(param_file), param_file)

    simulation = Simulation(space, space.star)
    simulation.run()