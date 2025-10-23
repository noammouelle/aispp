# contains basic functions needed to analyse example
# output files.
import h5py
import pandas as pd

def load_data(filename):
    with h5py.File(filename,"r") as file:
        states = file["states"][:]
        positions = file["positions"][:]
        velocities = file["velocities"][:]
        phase_shifts = file["phaseShifts"][:]
        phase_shifts_err = file["phaseShiftErrors"][:]
        interference_flag = file["interferingFlag"][:]

    df = pd.DataFrame({"states":states, "x":positions[:,0], "y":positions[:,1], "z":positions[:,2],
                    "vx":velocities[:,0], "vy":velocities[:,1], "vz":velocities[:,2],
                    "phase_shifts":phase_shifts, "phase_shifts_errors":phase_shifts_err,
                    "interference_flag":interference_flag})
    return df