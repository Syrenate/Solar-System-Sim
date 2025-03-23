import math
import numpy as np

class Vector2D:
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
        return np.arccos(np.clip(self.normalise().x, -1, 1))
    
    def normalise(self):
        mag = self.magnitude()
        return Vector2D((self.x / mag, self.y / mag))
    
class CelestialBody:
    """Class which represents an abstract celestial body in the solar system, i.e. a star or planet."""
    def __init__(self, reference_space, body_parameters, body_type):
        """Constructor method which sets the physical and visual parameters of the body."""
        self.name = body_parameters['name']
        self.mass = body_parameters['mass']
        self.colour = body_parameters['colour']
        self.type = body_type
        self.space = reference_space
        
        self.position = Vector2D((body_parameters['orbital_radius'], 0))
        # Initialise velocity based on Kepler approximation for a planets's stable orbit radius.
        init_velocity = math.sqrt(self.space.grav_constant * (self.space.star.mass) / body_parameters['orbital_radius']) if body_type == "planet" else 0
        self.velocity = Vector2D((0, init_velocity))
        self.acceleration = Vector2D((0,0))
        self.prev_acceleration = Vector2D((0,0))
        
        self.position_trail = []

    def beeman_update(self, is_final):
        """This method is called to update the position of the body, and subsequently apply the Beeman integration method to all bodies if all positions have been updated (is_final)."""
        # First predict the position of the body based on data at time steps (t) and (t - time_step).
        self.prev_position = self.position
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
        self.acceleration = self.get_acceleration(bodies_copy)
        self.velocity += self.acceleration.scale(self.space.time_step)
        self.position += self.velocity.scale(self.space.time_step)
        
    def direct_euler_update(self, bodies_copy):
        self.acceleration = self.get_acceleration(bodies_copy)
        self.position += self.velocity.scale(self.space.time_step)
        self.velocity += self.acceleration.scale(self.space.time_step)

    def get_acceleration(self, focus_objects):
        total_force = Vector2D((0,0))
        for focus in focus_objects:
            if self.name != focus.name:
                distance_sqr = (self.position.x - focus.position.x) ** 2 + (self.position.y - focus.position.y) ** 2
                total_force += (focus.position - self.position).normalise().scale((focus.mass) / distance_sqr)
        return total_force.scale(self.space.grav_constant)
    
class Planet(CelestialBody):
    """Class for planet to handle orbital periods and such."""
    def __init__(self, space, body_parameters, nasa_prediction):
        """Constructor method for planet, which defined the unique properties of orbit around a star."""
        super().__init__(space, body_parameters, "planet")

        self.orbital_period = 0
        self.nasa_orbital_period = nasa_prediction
        self.orbit_output = ""
        self.has_revolved = False

        self.angle_from_xaxis = 0

    def gen_orbit_output(self):
        """Method to generate a string that contains the essential orbit information of the planet."""
        self.orbit_output = f"     Planet: {self.name[0].upper() + self.name[1:] if len(self.name) > 1 else self.name.upper()}, Orbital Period: {round(self.orbital_period, 3)}Y"
        self.orbit_output += ((60 - len(self.orbit_output)) * ' ') + f"NASA Result: {self.nasa_orbital_period}Y\n"
