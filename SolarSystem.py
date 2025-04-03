import os
import json
import math
import time
import numpy as np

import matplotlib.path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Slider
from matplotlib.animation import FuncAnimation

from BodyClasses import CelestialBody, Planet

parameter_file_name = "simulation_parameters.json"
energy_file_name = "energy_level_data.txt"

class Space:
    """Class containing a complete physical model of the solar system."""
    def __init__(self, solar_parameters):
        """Constructor method in which the solar parameters are parsed to initialise the many-body system."""
        self.time_step = solar_parameters['timestep']               # Time step; the time interval to be simulated each frame.
        self.original_time_step = solar_parameters['timestep']      # Copy of the time_step parameter for use once the main value is changed.
        self.grav_constant = solar_parameters['grav_const']

        self.greatest_orbital_radii = 0    # Records greatest orbit to set axis limits in the simulaiton.
        self.elapsed_time = 0

        self.k_energy_level_history = []           # Record of previous energy levels for study.
        self.p_energy_level_history = []
        self.energy_level_history = []

        self.energy_history_sample_size = solar_parameters['experiment_2']['energy_sample_size']  # Desired maximum length of energy_level_history.
        self.energy_sample = solar_parameters['experiment_2']['energy_sample']                    # When recording the energy level, only sample every energy_sample frames.
        self.integration_type = solar_parameters['experiment_2']['integration_method']

        self.alignment_history = []      # Record occourances of planetary alignment.
        self.alignment_sample = solar_parameters['experiment_4']['planet_sample_size']        # How many planets should be considered for alignment.
        self.alignment_tolerance = solar_parameters['experiment_4']['alignment_tolerance'] * math.pi / 180    # Greatest angle between a planet and the mean angle to be considered alligned.
 
        self.planets = []       # Array of planets in the system.
        self.body_count = 0     # Total number of celestial bodies.
        for body in solar_parameters['bodies']:     # This loop scans all bodies specified in the json file, and creates instances to represent them.
            if body['type'] == "star": self.star = CelestialBody(self, body, "star", [(0,0), (0,0), (0,0)])    # self.star will reference the sun in our solar system.
            elif body['type'] == "planet": self.planets.append(Planet(self, body))
            self.body_count += 1 

            # Compare radius to the greatest radius recorded; replace if greater. This ensures all planets are on screen at 100% zoom.
            if body['orbital_radius'] > self.greatest_orbital_radii: self.greatest_orbital_radii = body['orbital_radius']

        for planet in self.planets:     # Once the greatest radius has been evaluated, set the trail sample of each planet. This should be proportional to the planet's orbital radius.
            planet.trail_length_limit = math.ceil(300 * (planet.position.x / self.greatest_orbital_radii) + 50)
  
        self.bodies = [self.star] + self.planets    # Array of all celestial bodies

    def update(self, iteration: int):
        """Method to simulate planetary motion in one time step."""
        self.determine_alignment()
        bodies_copy = self.bodies.copy()

        for i, body in enumerate(self.bodies):   # This loop updates all bodies in space, and evaluates if a planet has orbited or not.
            prev_pos = body.position
            if self.integration_type == 'beeman_integration': body.beeman_update(i == len(self.bodies) - 1)
            elif self.integration_type == 'euler_cromer': body.euler_cromer_update(bodies_copy)
            elif self.integration_type == 'direct_euler': body.direct_euler_update(bodies_copy)

            if body.type == "planet":
                # If not the first frame, the planet has not revolved before, the planet is in the +ve x-axis region, and it crosses the axis boundary; it has revolved. Record orbital period.
                if iteration != 0 and not body.has_revolved and body.position.x > 0 and (body.position.y > 0 and prev_pos.y <= 0):
                    body.orbital_period = self.elapsed_time    # Orbital period is the time since the simulation started.
                    body.has_revolved = True    
                    body.generate_orbit_output()

            if body.type != "star":     # Every non-star body should have a trail (since the star should be centered).
                if iteration % body.trail_sample == 0: body.position_trail.append((body.position - self.star.position).get())
                if len(body.position_trail) > body.trail_length_limit: body.position_trail = body.position_trail[1:]

        if len(self.energy_level_history) > self.energy_history_sample_size: 
                self.energy_level_history = self.energy_level_history[1:]
        self.elapsed_time += self.time_step
    
    def determine_alignment(self):
        """Method which checks planetary alignment"""
        alignments = []
        for planet in self.planets[:self.alignment_sample]:    # First generate the angle from 0 to 2pi with the +ve x-axis for each planet.
            alignments.append(planet.position.angle_between_xaxis(self.star.position))

        mean_angle = np.mean(alignments)
        alligned = True     
        # Now loop through each plannet to check if they are close enough to the mean angle. Is so, they are 'alligned'.
        for i, planet in enumerate(self.planets[:self.alignment_sample]):
            if alligned:
                angle = alignments[i]
                angle_min, angle_max = min([angle, mean_angle]), max([angle, mean_angle])
                # Find the minimum angle between the planet and the mean angle, and evaluate if this is within the threshold for alignment.
                min_angle_delta = min([angle_max - angle_min, 2 * math.pi - (angle_max - angle_min)])
                if min_angle_delta > self.alignment_tolerance: alligned = False

        # If planets allign this frame, only append to the list of occourances of alignment if we're certain the current alignment has not been recorded.
        if (len(self.alignment_history) == 0 and alligned) or (self.planets[0].has_revolved and alligned and np.abs(self.alignment_history[-1][0] - self.elapsed_time) > self.planets[0].orbital_period): 
            self.alignment_history.append((float(self.elapsed_time), min_angle_delta * (180 / math.pi)))

    def get_energy_levels(self):
        """Method that calculates and returns the kinetic & potential energy levels of all bodies."""
        k_energy, p_energy = 0, 0
        for body in self.bodies:
            k_energy += body.mass * (body.velocity.magnitude() ** 2)
            p_energy -= 0 if body.type == "star" else body.mass / (body.position - self.star.position).magnitude()
        return k_energy / 2, p_energy * (self.star.mass * self.grav_constant)
    
class Simulation:
    """Class to represent the matplotlib simulation space, handles UI elements and simulation rendering."""
    def __init__(self, reference_space: Space, parameters):
        """Constructor method for the simulation; defines the simulation and UI axes in matplotlib, according to a reference Space instance."""
        self.paused = False
        self.default_os = parameters['default_os']              # The operating system specified; used for clearing the console.
        self.iteration_limit = parameters['num_iterations']     # Maximum iterations for the simulation. Terminates once reached.

        self.elapsed_frames = 0      
        # To estimate framerate, record the time taken for recent frames and take an average.   
        self.framerate_history = []          # Store recent frame data.
        self.framerate_history_length = 100  # Maximum sample size for estimation.
        
        colour = 'white'
        plt.rcParams['text.color'] = colour
        plt.rcParams['axes.labelcolor'] = colour
        plt.rcParams['xtick.color'] = colour
        plt.rcParams['ytick.color'] = colour
        self.body_sizes = {"star": 0.2, "planet": 0.1, "moon": 0.05}

        self.space = reference_space        # Keep reference to the instance of Space we're simulating.
        self.space_reference_object = reference_space.star
        self.sim_fig, self.sim_ax = plt.figure(1), self.generate_simulation_axes()   # Initialise the simulation axes.

        self.sim_fig.set_facecolor((0.12,0.12,0.12))
        self.sim_ax.patch.set_edgecolor('white')
        self.sim_ax.patch.set_linewidth(2)
        self.sim_ax.set_facecolor((0.02,0.02,0.02))
        plt.title("Many-Body Simulation of our Solar System", color='white')

        def update_simulation_limit(val):   # Method to rescale the simulation boundary to zoom in/out.
            scale = self.space.greatest_orbital_radii / val
            self.sim_ax.set_xlim(-scale, scale)
            self.sim_ax.set_ylim(-scale, scale)

        def update_time_scale(val):     # Method which updates the space's time scale.
            self.space.time_step = self.space.original_time_step * val

        self.space.current_energy_level = sum(self.space.get_energy_levels())

        self.sim_ax.axis("scaled")
        update_simulation_limit(5)
        self.sim_fig.subplots_adjust(left=0.3, bottom = 0.3)

        axzoom = self.sim_fig.add_axes([0.35, 0.15, 0.55, 0.03], xscale='log')
        self.zoom_slider = Slider(      # Interactive slider to zoom in/out of the simulation space.
            ax=axzoom,
            label='Zoom scale',
            valmin=0.9,
            valmax=10,
            valinit=5,
        )

        axtime = self.sim_fig.add_axes([0.35, 0.05, 0.55, 0.03], xscale='log')
        self.time_slider = Slider(      # Interactive time scale slider for the space instance.
            ax=axtime,
            label='Time scale',
            valmin=0.1,
            valmax=15,
            valinit=1,
        )
        def on_click(e): self.paused = True
        def on_release(e): self.paused = False
        self.sim_fig.canvas.mpl_connect('button_press_event', on_click)
        self.sim_fig.canvas.mpl_connect('button_release_event', on_release)

        self.zoom_slider.on_changed(update_simulation_limit)
        self.time_slider.on_changed(update_time_scale)

        self.start_update_time = time.time()

    def update(self, iteration):
        """Method that handles iterating to the next frame of the simulation."""
        if self.paused: return []

        self.elapsed_frames += 1

        self.space.update(iteration)    # Update the reference space immediately.
        for i in range(len(self.space.bodies)):     # Modify the patches of planets and their trails based on this update.
            body = self.space.bodies[i]
            self.patches[i].center = (body.position - self.space_reference_object.position).get()
            if body.type == "planet": self.patches[i + self.space.body_count - 1].set_path(matplotlib.path.Path(body.position_trail))

        # Calculation for framerate
        self.framerate_history.append(time.time() - self.start_update_time)
        if len(self.framerate_history) > self.framerate_history_length: self.framerate_history = self.framerate_history[1:]

        k_energy, p_energy = self.space.get_energy_levels()
        if iteration % self.space.energy_sample == 0:
            self.space.k_energy_level_history.append(k_energy)
            self.space.p_energy_level_history.append(p_energy)
            self.space.energy_level_history.append(k_energy + p_energy)
            
            with open(energy_file_name, 'a') as file: 
                time_str = f"Frame {iteration} ({round(self.space.elapsed_time,2)}Y)"
                energy_str = f"Energy levels: [Kinetic: {round(k_energy,3)}, Potential: {round(p_energy,3)}]."
                file.write(time_str + (' '*(25 - len(time_str))) + (energy_str + (' '*(60 - len(energy_str)))) + f"Total energy: {round(k_energy + p_energy,3)}.\n")
                            
        # Generate a console output of simulation data every 30 frames.
        if iteration % 30 == 0 or iteration == self.iteration_limit - 1:
            if self.default_os == "windows": os.system('cls')
            elif self.default_os == "linux": os.system('clear')
            print(self.generate_output(k_energy, p_energy))

        self.start_update_time = time.time()
        return self.patches         # FuncAnimation requires returning all artists which should be rendered.

    def run(self):
        """Method called to initialise the simulation."""
        self.simulation_start_time = time.time()
        self.anim = FuncAnimation(self.sim_fig, self.update, frames=self.iteration_limit, repeat=False, interval=1/30000, blit=True)
        plt.show()

    def generate_simulation_axes(self):
        """Method to create patches for each simulation element; all physical objects and their trails if applicable."""
        ax = plt.axes()
        self.patches = []
        for body in [self.space.star] + self.space.planets:     # Render each celestial body
            # Planets further from the star should have greater radius for better visibility.
            #radius_scalar = 0.1 + 0.3 * ((body.position.x / self.space.greatest_orbital_radii) ** 4) if body.type == "planet" else 0.1 
            self.patches.append(mpatches.Circle(body.position.get(), self.body_sizes[body.type], color = body.colour, animated = True))
            ax.add_patch(self.patches[-1])

        for planet in self.space.planets:     # Render the trail of each planet
            planet.position_trail.append(planet.position.get())
            path = matplotlib.path.Path(planet.position_trail, codes = [np.uint8(1)], closed=True)
            self.patches.append(mpatches.PathPatch(path, color = planet.colour, animated = True, fill = False))
            ax.add_patch(self.patches[-1])

        return ax

    def generate_output(self, k_energy, p_energy):
        """Method which generates a console output of simulation and experiment data."""
        output = ""

        elapsed_real_time = (time.time() - self.simulation_start_time) / 60
        frame_rate = 1/np.mean(self.framerate_history)

        simulation_output =  f"     Elapsed simulation time: {round(self.space.elapsed_time, 3)} Years  [Real time: {math.floor(elapsed_real_time)} min {round((elapsed_real_time - math.floor(elapsed_real_time)) * 60)}s]\n"
        simulation_output += f"     Current frame: {self.elapsed_frames}/{self.iteration_limit}"
        simulation_output += f" ({round(100 * self.elapsed_frames / self.iteration_limit, 2)}% elapsed, framerate: {round(frame_rate, 2)}fps)\n"
        simulation_output += f"     Simulation speed: {round(self.space.time_step * frame_rate, 2)} years per simulation second."
        output += "Simulation data:\n" + simulation_output + "\n\n"

        # Experiment 1
        planet_output = ""
        for planet in self.space.planets: planet_output += planet.orbit_output
        if len(planet_output) > 0: output += "Orbital Periods (Experiment 1):\n" + planet_output + "\n"

        # Experiment 2
        energy_output =  f"     Kinetic Energy: {round(k_energy, 3)}, Potential Energy: {round(p_energy, 3)}\n"
        energy_output += f"     Total Energy: {round(k_energy + p_energy, 3)}. "
        energy_output += f"Recent standard deviation [sample size: {len(self.space.energy_level_history)}]: {round(math.sqrt(np.var(self.space.energy_level_history)), 3)}"
        output += "Energy Levels (Experiment 2):\n" + energy_output + "\n\n"

        # Experiment 4
        alignmnet_output = f"Alignment occourances (Experiment 4): \n"
        alignmnet_output += f"     - Current study: {self.space.alignment_sample} planets align within {round(self.space.alignment_tolerance * 180 / math.pi)} degrees of their mean.\n"
        for data in self.space.alignment_history: 
            alignmnet_output += f"         > At {round(data[0], 2)}Y, planets align within {round(data[1], 2)} degrees of their mean.\n"
        output += alignmnet_output

        return output

if __name__ == "__main__":
    with open(energy_file_name, 'w') as file: file.write("")
    with open(parameter_file_name) as file: 
        parameters = json.load(file)
        space = Space(parameters)

    simulation = Simulation(space, parameters)
    simulation.run()