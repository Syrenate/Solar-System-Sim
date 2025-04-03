import math
import numpy as np

class Vector2D:
    """Class to represent and handle arithmetic of 2D vectors, used as acceleration vectors for example."""
    def __init__(self, pos):
        self.set(pos)

    def set(self, pos):
        self.x, self.y = pos[0], pos[1]

    def get(self):
        return (self.x, self.y)

    def __add__(self, other):
        return Vector2D((self.x + other.x, self.y + other.y))
    
    def __sub__(self, other):
        return Vector2D((self.x - other.x, self.y - other.y))
    
    def scale(self, scalar):
        return Vector2D((self.x * scalar, self.y * scalar))
    
    def magnitude(self):
        return math.sqrt(self.x ** 2 + self.y ** 2)
    
    def dot(self, pos):
        return self.x * pos.x + self.y * pos.y

    def angle_between_xaxis(self):
        """Method to return the angle between the vector and the positive x-axis, specifically the vector (1,0)."""
        angle = np.arccos(np.clip(self.normalise().x, -1, 1))
        if self.y < 0: return 2 * math.pi - angle
        else: return angle
    
    def normalise(self): 
        return Vector2D((self.x, self.y)).scale(1/self.magnitude())
    
class Object:
    """Class representing a physical object with momentum in some space."""
    def __init__(self, space, parameters, init):
        self.name = parameters['name']
        self.mass = parameters['mass']
        self.space = space

        self.position = Vector2D(init[0])
        self.velocity = Vector2D(init[1])
        self.acceleration = Vector2D(init[2])
        self.prev_acceleration = Vector2D((0,0))   # For use in Beeman integration.

    def beeman_update(self, is_final):
        """Beeman integration method to update positions, and subsequently update velocity of all bodies if all positions have been updated (is_final)."""
        # First predict the position of the body based on data at time steps (t) and (t - time_step).
        self.position += self.velocity.scale(self.space.time_step) + (self.acceleration.scale(4) - self.prev_acceleration).scale(((self.space.time_step ** 2) / 6))

        if is_final:    # If all positions have been updatad, evaluate new velocity of each body.
            for i in range(len(self.space.bodies)):
                body = self.space.bodies[i]

                # Compute acceleration based on new positions of bodies. Only celestial bodies will impact a bodies's acceleration meaningfully.  
                new_acceleration = body.get_acceleration(body.space.bodies)
                body.velocity += (new_acceleration.scale(2) + body.acceleration.scale(5) - body.prev_acceleration).scale(body.space.time_step / 6)

                body.prev_acceleration = body.acceleration
                body.acceleration = new_acceleration

    def euler_cromer_update(self, bodies_copy):
        """Method to implement Euler-Cromer integration for updating object positions."""
        self.acceleration = self.get_acceleration(bodies_copy)
        self.velocity += self.acceleration.scale(self.space.time_step)
        self.position += self.velocity.scale(self.space.time_step)
        
    def direct_euler_update(self, bodies_copy):
        """Method to implement direct-Euler integration for updating object positions."""
        self.acceleration = self.get_acceleration(bodies_copy)
        self.position += self.velocity.scale(self.space.time_step)
        self.velocity += self.acceleration.scale(self.space.time_step)

    def get_acceleration(self, focus_objects):
        """Method to get the acceleration of self under the influence of a list of focus_bodies."""
        acceleration = Vector2D((0,0))
        for focus in focus_objects:
            if self.name != focus.name:
                distance_sqr = (self.position.x - focus.position.x) ** 2 + (self.position.y - focus.position.y) ** 2
                acceleration += (focus.position - self.position).normalise().scale((focus.mass) / distance_sqr)
        return acceleration.scale(self.space.grav_constant)

class CelestialBody(Object):
    """Class which represents an abstract celestial body in the solar system, i.e. a star or planet."""
    def __init__(self, reference_space, body_parameters, body_type, init):
        """Constructor method which sets the physical and visual parameters of the body."""
        super().__init__(reference_space, body_parameters, init)
        self.colour = body_parameters['colour']
        self.type = body_type
        
        self.trail_sample = 2           # Sample position every this many frames.
        self.position_trail = []        # Record of recent positions to create a path between them, creating a trail of the body.
    
class Planet(CelestialBody):
    """Class for planet to handle orbital periods and such."""
    def __init__(self, reference_space, body_parameters):
        """Constructor method for planet, which defines the properties for orbiting a star."""

        # Initialise position and velocity based on Kepler approximation for a planets's stable orbit radius.
        init_pos, init_velocity = body_parameters['orbital_radius'], math.sqrt(reference_space.grav_constant * (reference_space.star.mass) / body_parameters['orbital_radius'])
        super().__init__(reference_space, body_parameters, "planet", [(init_pos, 0), (0,init_velocity), (0, 0)])
    
        self.orbital_period = 0
        self.is_study = body_parameters['is_study']
        self.nasa_orbital_period = body_parameters['nasa_orbit_period'] if self.is_study else 0
        self.orbit_output = ""
        self.has_revolved = False

        self.angle_from_xaxis = 0

    def generate_orbit_output(self):
        """Method to generate a string that contains the essential orbit information of the planet."""
        self.orbit_output = f"     Planet: {self.name[0].upper() + self.name[1:] if len(self.name) > 1 else self.name.upper()}, Orbital Period: {round(self.orbital_period, 3)}Y"
        self.orbit_output += ((60 - len(self.orbit_output)) * ' ') + (f"NASA Result: {self.nasa_orbital_period}Y [delta: {round(self.nasa_orbital_period - self.orbital_period, 3)}]\n" if self.is_study else "\n")