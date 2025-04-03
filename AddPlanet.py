import json
import time

if __name__ == "__main__":
    run_interface = True
    while run_interface:
        print("Welcome to the console interface for simulation modification. Please select an option:")
        print("     1. Create a new planet")
        print("     2. Modify an existing planet")
        print("     3. Remove a planet")
        print("     4. Exit the interface.\n")
        user_choice = input("Please select from options 1-4: "); time.sleep(0.5)

        if user_choice[0] in ['2','3']: 
            print("Current bodies: ")
            with open("simulation_parameters.json", 'r') as file:
                bodies = json.load(file)['bodies']
                for body in bodies:
                    print(f"    -{body['type']}, name: {body['name']}.")

        if user_choice[0] == '1':
            print("Creating planet. Please ensure units are SI relative to Earth (for example, mass should be the ratio of YOUR planet's mass to EARTH'S.).")
            planet = {}
            planet["name"] = input("Planet name: "); time.sleep(0.3)
            planet["type"] = "planet"
            planet["mass"] = float(input("Planet's mass: ")); time.sleep(0.3)
            planet["orbital_radius"] = float(input("Orbital radius: ")); time.sleep(0.3)
            planet["is_study"] = False
            planet["colour"] = input("Planet's colour: "); time.sleep(0.5)

            with open("simulation_parameters.json", 'r') as file:
                parameters = json.load(file)
                bodies = parameters['bodies']
                for i, body in enumerate(bodies):
                    if body['orbital_radius'] > planet['orbital_radius']:
                        bodies = bodies[:i] + [planet] + bodies[i:]; break
                parameters['bodies'] = bodies
            with open("simulation_parameters.json", 'w') as file:
                json.dump(parameters, file)
            
        if user_choice[0] == '2':
            name = input("Modifying planet. Enter the planet's name: ")
            time.sleep(0.3)

            found = False
            with open("simulation_parameters.json", 'r') as file:
                parameters = json.load(file)
                bodies = parameters['bodies']
                for i, body in enumerate(bodies):
                    if body['name'] == name:
                        planet = body
                        planet_index = i
                        found = True
                        break

            if found == False:
                print("Planet not found!")
                break
            else: 
                choice = input(f"Enter the property you wish to modify {planet.keys()}: "); time.sleep(0.3)
                value = input(f"Enter the value to change to: [current: {planet[choice] }]"); time.sleep(0.3)

                if not choice in ["name", "colour", "type"]: value = float(value)
                planet[choice] = value
                parameters['bodies'] = bodies[:planet_index] + [planet] + bodies[planet_index:]

                with open("simulation_parameters.json", 'w') as file:
                    json.dump(parameters, file)

                print("Modified!")


        if user_choice[0] == '3':
            name = input("Removing planet. Enter the planet's name: ")

            found = False
            with open("simulation_parameters.json", 'r') as file:
                parameters = json.load(file)
                bodies = parameters['bodies']
                for i, body in enumerate(bodies):
                    if body['name'] == name:
                        bodies = bodies[:i] + bodies[i+1:]
                        found = True
                        break
                parameters['bodies'] = bodies

            if found == False:
                print("Planet not found!")
            else:
                with open("simulation_parameters.json", 'w') as file:
                    json.dump(parameters, file)
                time.sleep(0.5)
                print("Removed!")
        elif user_choice[0] == '4': run_interface = False
        print("\n")