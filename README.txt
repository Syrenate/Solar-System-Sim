Hello! This repo contains all necessary components for running and experimenting with a Solar System simulation.

Requirements: Windows 10 System, Matplitlib, Numpy.
    - If using Linux, modify the value of 'default_os' to 'linux' in simulation_parameters.json.

Instructions for use:
    - Running SolarSystem.py with the most recent Python version will render a complete simulation with most experiment data available without modification.
    - Experiment 2 investigates different integration techniques; this may be modified in simulation_parameters.json under 'integration_method'. 
        - Options: 'beeman_integration', 'euler_cromer', or 'direct_euler'.
    - Most experimental data is output to the console, including some simulation data.
    - Interactive matplotlib sliders allow the real-time adjustment of certain parameters.

    - To add, modify or remove a planet, run AddPlanet.py and follow console instructions.