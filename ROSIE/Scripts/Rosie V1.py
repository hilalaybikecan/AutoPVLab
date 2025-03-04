import time
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

current_dir = Path(__file__)
# Find folder "rosie", which has to contain "uarm"
project_dir = str([p for p in current_dir.parents if p.parts[-1] == "rosie"][0])
sys.path.append(project_dir)

# sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
from uarm.wrapper import SwiftAPI

import serial


"""
api test: move
"""


def convert_ml_to_steps(ml):
    return int(np.round(5625*ml)) #5625 is for 0.9 ml becasue syringe is not allowed to mowe for 1 ml yet.

def cartesian_to_polar_coordinates(arr, extension=35):
    """ Converts cartesian to polar coordinates.
    Parameters
    -----------
    arr : np.array
    extension : float
        Distance between pipet holder and default hand coordinate.
    """
    x = arr[:, 0]  # 0-th column
    y = arr[:, 1]  # 1-st column

    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)/np.pi*180 #Element-wise arc tangent of x1/x2 choosing the quadrant correctly.
    # The quadrant (i.e., branch) is chosen so that arctan2(x1, x2) is the signed angle in radians between the ray
    # ending at the origin and passing through the point (1,0), and the ray ending at the origin and passing through
    # the point (x2, x1). (Note the role reversal: the “y-coordinate” is the first function parameter,
    # the “x-coordinate” is the second.)

    # Make sure, coordinates are legal / valid
    assert np.all(r-extension > 95), "Too close to Rosie"
    # assert np.all(phi >= 0) and np.all(phi <= 180), "Angles out of range"

    return np.round(np.array([r-extension, 90 + phi]).T, decimals=4)


def polar_to_cartesian_coordinates(arr, extension=35):
    """ Converts cartesian to polar coordinates.
    Parameters
    -----------
    arr : np.array
    extension : float
        Distance between pipet holder and default hand coordinate.
    """
    r = arr[:, 0] + extension # 0-th column
    phi = arr[:, 1]  # 1-st column

    # The quadrant (i.e., branch) is chosen so that arctan2(x1, x2) is the signed angle in radians between the ray
    # ending at the origin and passing through the point (1,0), and the ray ending at the origin and passing through
    # the point (x2, x1). (Note the role reversal: the “y-coordinate” is the first function parameter,
    # the “x-coordinate” is the second.)
    x = r*np.cos((phi-90)/180*np.pi)
    y = r*np.sin((phi-90)/180*np.pi)

    # Make sure, coordinates are legal / valid

    return np.array([x, y]).T


def generate_grid(x_0=0, y_0=0, x_N=8, y_N=12, Delta_x=10, Delta_y=10):
    """
    Parameters
    ----------
    x_0, y_0 : float
        Upper left corner of grid (not from Rosie's perspective)
    x_N, y_N : int
        Number of rows and columns
    Delta_x, Delta_y : float
        Distance between rows or columns
    """
    grid = [np.array([x, y]) + np.array([x_0, y_0])
            for x in np.arange(x_N)*Delta_x for y in np.arange(y_N)*Delta_y]
    return np.array(grid)


def get_pipet(r, phi, z_lower=50, z_upper=80, speed=10_000):
    """ Move to coordinates (x,y), lower arm to z, get pipet, move up again.
    Parameter
    ---------
    x : float
        X-coordinate
    y : float
    """
    # Move into position of pipet
    swift.set_polar(stretch=r, rotation=phi, height=z_upper, speed=speed*50, wait=True, timeout=None)
    swift.flush_cmd(wait_stop=True)

    # Go down
    # swift.set_position(z=z, speed=speed, wait=True)
    swift.set_polar(stretch=r, rotation=phi, height=z_lower, speed=speed, wait=True, timeout=None)
    time.sleep(0.25)

    # Go up again
    swift.set_polar(stretch=r, rotation=phi, height=z_upper, speed=speed, wait=True, timeout=None)
    swift.flush_cmd(wait_stop=True)


def handle_solution(syringe_val=0.05, with_syringe=True, pause=4):
    if with_syringe:
        arduino_msg = f"{int(convert_ml_to_steps(syringe_val))}"  # format integer as a string
        arduino_msg = arduino_msg.encode("ascii")  # encode string in utf-8 (that is no special characters)

        # Send message to Arduino
        arduinoData.write(arduino_msg)
        arduinoData.flush()
        arduinoData.reset_input_buffer()
        arduinoData.reset_output_buffer()

        # time.sleep(5)
        time.sleep(abs(syringe_val)*10*5 + pause)

        # arduino_has_finished = False
        # while not arduino_has_finished:
        #     number_of_steps_inputted = arduinoData.readline()
        #     arduino_has_finished = (number_of_steps_inputted == b"done\r\n")
        #     print(arduino_has_finished)
        #     print("Waiting")
        #     time.sleep(0.2)


        # Tell us how much the arduino moved
        # number_of_steps_inputted = arduinoData.readline()
        # print(number_of_steps_inputted)


def go_mix_pipet(r, phi, z_lower=25, z_upper=50, speed=10_000, syringe_val=0.05, with_syringe=True, pause=4):
    """ Move to coordinates (x,y), lower arm to z, empty pipet, move up again.
    Parameter
    ---------
    x : float
        X-coordinate
    y : float
    syringe_val : float, in ml

    """
    # Move into position of pipet
    swift.set_polar(stretch=r, rotation=phi, height=z_upper, speed=speed, wait=True, timeout=None)
    swift.flush_cmd(wait_stop=True)

    # Go down
    # swift.set_position(z=z, speed=speed, wait=True)
    swift.set_polar(stretch=r, rotation=phi, height=z_lower, speed=speed, wait=True, timeout=None)

    time.sleep(0.5)

    # Communication with Arduino
    handle_solution(syringe_val=syringe_val, with_syringe=with_syringe, pause=pause)

    # Go up again
    swift.set_polar(stretch=r, rotation=phi, height=z_upper, speed=speed, wait=True, timeout=None)
    swift.flush_cmd(wait_stop=True)



def put_into_trash(x, y, z=150, speed=1e6):
    # Move into position of pipet
    swift.set_position(x=x, y=y, z=z, speed=speed, wait=True, timeout=None)

    time.sleep(2)

    # Drop pipet
    swift.set_wrist(180, wait=True, timeout=None, speed=speed)
    time.sleep(0.25)
    swift.set_wrist(90, wait=True, timeout=None, speed=speed)
    swift.flush_cmd(wait_stop=True)


def distribute_load_to_pipets(n_pipets, total_sol, comb_sol, mask, round_to=3):
    # Initialise return lists
    coord_per_pipet = [[] for i in range(n_pipets)]
    syringe_val_per_well = [[] for i in range(n_pipets)]
    if total_sol/0.2 - np.floor(total_sol/0.2) > 0:
        amount_in_pipet = (
            [0.2 for i in range(n_pipets-1)]
            + [np.round((total_sol/0.2 - np.floor(total_sol/0.2))*0.2, 10)]
        )
    else:
        amount_in_pipet = [0.2 for i in range(n_pipets)]

    amount_in_pipet_return = amount_in_pipet.copy()
    comb_sol_tmp = comb_sol.copy()

    k_pipet = 0  # iterates over pipets
    j_well = 0  # well
    while j_well < len(comb_sol):
        coord_per_pipet[k_pipet].append(mask[j_well])

        if np.round((amount_in_pipet[k_pipet] - comb_sol_tmp[j_well]), round_to) > 0:  # we don't need to use several pipets
            amount_in_pipet[k_pipet] = np.round(amount_in_pipet[k_pipet] - comb_sol_tmp[j_well], round_to)
            syringe_val_per_well[k_pipet].append(np.round(comb_sol_tmp[j_well], round_to))
            comb_sol_tmp[j_well] = 0.0
        else:  # Not enough in the current pipet, so we need the next pipet to come back to this well
            syringe_val_per_well[k_pipet].append(np.round(amount_in_pipet[k_pipet], round_to))  # empty pipet
            comb_sol_tmp[j_well] = np.round(comb_sol_tmp[j_well] - amount_in_pipet[k_pipet], round_to)  #
            amount_in_pipet[k_pipet] = 0.0
            k_pipet += 1  # go to next pipet
            assert k_pipet <= n_pipets, "Not enough pipets"
        if comb_sol_tmp[j_well] == 0:
            j_well += 1

    assert np.all(np.round(amount_in_pipet, 10) == 0), "Logic error in distribution of pipets!"
    return coord_per_pipet, syringe_val_per_well, amount_in_pipet_return

# Whether to try send commands to the Arduino controlling the syringe
with_syringe = True# False



# Initialises class instant swift that lets us control Rosie
swift = SwiftAPI(filters={'hwid': 'USB VID:PID=2341:0042'})
# Initialise comunication instant for Arduino controlling the syringe
if with_syringe:
    arduinoData = serial.Serial("COM4", 9600, write_timeout = 0)

    while (arduinoData.inWaiting() == 0):
        pass
    # Is the arduino talking to us?
    test_string = arduinoData.readline()
    print(test_string)


# Wait until Rosie is ready
swift.waiting_ready(timeout=3)

device_info = swift.get_device_info()
print(device_info)
firmware_version = device_info['firmware_version']
if firmware_version and not firmware_version.startswith(('0.', '1.', '2.', '3.')):
    swift.set_speed_factor(0.0005)

swift.set_mode(0)  # general mode

#print(swift.get_position())
# quit()
swift.reset(wait=True, speed=50000, x=190, y=0, z=160)
swift.flush_cmd(wait_stop=True)
#time.sleep(3)








# Main code starts here

# Define grids
pipet_grid = cartesian_to_polar_coordinates(
    generate_grid(x_0=35, y_0=-242, x_N=6, y_N=11, Delta_x=-8.9, Delta_y=-8.9),
    extension=35)

solution_grid = cartesian_to_polar_coordinates(
    generate_grid(x_0=78.5, y_0=-284.5, x_N=2, y_N=4, Delta_x=25, Delta_y=25),
    extension=35)

# x = kisa kenar, y = uzun kenar

x_N_well, y_N_well = 5, 6
well_plate_grid = cartesian_to_polar_coordinates(
    generate_grid(x_0=150, y_0=-246.3, x_N=x_N_well, y_N=y_N_well, Delta_x=8.9, Delta_y=8.9),
    extension=35)

print("Pipet coordinates")
print(pipet_grid)

# Invent a list of recepies
np.random.seed(42)
# combinations = (
#     (np.random.uniform(0, 1, (N_well_plate, N_sol)) < 0.7)
#     * np.round(np.random.uniform(0, 0.2, (N_well_plate, N_sol)), 2)
# )

# Read in csv file
df = pd.read_csv(r"C:\Users\can\Desktop\ROSIE\com_for_Rosie.csv", index_col=0)
df = df.dropna(axis=1, how='all', inplace=False)  # Delete columns (axis=1) that are empty / contain NaN values
#exclude = []  # ['CsCl', 'PbCl2']
#df = df.loc[:, [col for col in df.columns if not col in exclude]]
# df.iloc[:, 2:3]
combinations = df.to_numpy()[0:x_N_well*y_N_well, :]


print("Combinations to try")
print(combinations)
readable_coords = np.array([[chr(65+letter) + str(1+number)] for letter in np.arange(x_N_well)
                            for number in np.arange(y_N_well)])

#assert combinations.shape[0] <= well_plate_grid.shape[0], "Too many combinations specified in csv."
#assert combinations.shape[1] <= solution_grid.shape[0], "Too many solutions speficied in csv."
# N_well_plate = combinations.shape[0]  # well_plate_grid.shape[0]
N_sol = combinations.shape[1]  # solution_grid.shape[0]


# df = pd.DataFrame(
#     combinations,
#     columns=[f"Sol {i}" for i in range(N_sol)]
# )

for i_sol in range(0, N_sol):
    mask = np.where(combinations[:, i_sol] > 0)[0]
    comb_sol = combinations[mask, i_sol]

    total_sol = np.round(np.sum(combinations[:, i_sol]), 8)  # total amount of solution required
    n_pipets = int(np.ceil(total_sol/0.2))

    if total_sol == 0:
        continue  # skip this solution, as none is distributed

    coord_per_pipet, syringe_val_per_well, syringe_val_of_pipet = (
        distribute_load_to_pipets(n_pipets, total_sol, comb_sol, mask, round_to=3)
    )
    print(coord_per_pipet, syringe_val_per_well, syringe_val_of_pipet)

    # Get the pipet
    get_pipet(r=pipet_grid[i_sol, 0], phi=pipet_grid[i_sol, 1],
              z_lower=82, z_upper=160, speed=5e4)
    time.sleep(10)

    print("")  # for new line
    print(f"Start distributing Solution {i_sol} using {n_pipets} pipets.")
    for i_pipet in range(n_pipets):
        print(f"Load {syringe_val_of_pipet[i_pipet]} ml of Solution {i_sol} and distribute")
        for el in np.char.add(
                np.char.add(np.array(syringe_val_per_well[i_pipet], dtype=str).reshape(-1, 1), " ml to "),
                 readable_coords[coord_per_pipet[i_pipet]],
        ):
            print(el[0])
        print(" ")

        # Get solution
        step_to_red_error = 0.004  # in ml
        # print(convert_ml_to_steps(step_to_red_error))
        # 22 steps of pulling
        go_mix_pipet(r=solution_grid[i_sol, 0], phi=solution_grid[i_sol, 1],
                     z_lower=150, z_upper=150, speed=1e5,
                     syringe_val=step_to_red_error,
                     with_syringe=with_syringe, pause=1)
        # Pull solution
        go_mix_pipet(r=solution_grid[i_sol, 0], phi=solution_grid[i_sol, 1],
                     z_lower=83, z_upper=83, speed=1e5,
                     syringe_val=syringe_val_of_pipet[i_pipet],
                     with_syringe=with_syringe)
        # 22 steps of pushing
        go_mix_pipet(r=solution_grid[i_sol, 0], phi=solution_grid[i_sol, 1],
                     z_lower=83, z_upper=83, speed=1e5,
                     syringe_val=-step_to_red_error,
                     with_syringe=with_syringe)  # pause=1 is deleted here
        time.sleep(3)
        go_mix_pipet(r=solution_grid[i_sol, 0], phi=solution_grid[i_sol, 1],
                     z_lower=152, z_upper=152, speed=1e5,
                     syringe_val=-0,
                     with_syringe=False)

        # Distribute solution
        for i_coord in range(len(coord_per_pipet[i_pipet])):
            # TODO something else
            r, phi = well_plate_grid[coord_per_pipet[i_pipet][i_coord]]

            go_mix_pipet(r=r, phi=phi, z_lower=85, z_upper=130,  speed=1e5,
                         syringe_val=-syringe_val_per_well[i_pipet][i_coord], with_syringe=with_syringe)

    put_into_trash(x=150, y=0, speed=5e5)
    swift.reset(wait=True, speed=1e5)

# Sometimes I encountered the problem that it would skip one of the movements.
# This is almost always due to illegal coordinates being handed to Rosie.

swift.reset(wait=True, speed=1e5)

swift.flush_cmd()

time.sleep(3)
swift.disconnect()

# LAST UPDATE July 10, 2024  Hilal




