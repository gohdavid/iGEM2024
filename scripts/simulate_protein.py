#!/usr/bin/env python3

from openmm import app as app
from openmm import *
import mdtraj as md
import numpy as np
import time
import argparse
import os
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Process a PDB file and output results to a specified folder.')
    parser.add_argument('-i', '--input', required=True, help='Input PDB file name.')
    parser.add_argument('-o', '--output', required=True, help='Output folder path.')
    parser.add_argument('-t', '--time', required=True, help='Simulation time in ns.')
    parser.add_argument('-w', '--water', required=True, help='0 for implicit, 1 for explicit')
    
    args = parser.parse_args()

    input_file = Path(args.input)
    output_folder = Path(args.output)
    simulation_time = np.float64(args.time)
    water = bool(args.water)

    # Check if input file exists
    if not input_file.is_file():
        raise FileNotFoundError(f"The file {input_file} does not exist!")

    # Ensure the output directory exists, create if it does not
    if not output_folder.exists():
        output_folder.mkdir(parents=True)
        print(f"Created output directory: {output_folder}")

    # Load the PDB structure
    pdb = app.PDBFile(input_file.as_posix())
    modeller = app.Modeller(pdb.topology, pdb.positions)

    if water:
        # Add solvent with a specified ionic strength
        forcefield = app.ForceField('amber14/protein.ff14SB.xml','amber14/tip3p.xml')
        # Set up the system
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.9*unit.nanometers, constraints=app.HBonds,
                                        hydrogenMass=1.5*unit.amu, ewaldErrorTolerance=0.0005) # Use HMR
        # Add pressure
        barostat = MonteCarloBarostat(1*unit.bar,310.15*unit.kelvin)
        barostat_id = system.addForce(barostat)
    else:
        # Use implicit water
        modeller.deleteWater()
        ion_residues = ['NA', 'CL']
        # Find all ion residues to delete
        residues_to_delete = [residue for residue in modeller.topology.residues() if residue.name.upper() in ion_residues]
        # Remove the ion residues
        if residues_to_delete:
            modeller.delete(residues_to_delete)
        # Add solvent with a specified ionic strength
        forcefield = app.ForceField('amber14/protein.ff14SB.xml', 'implicit/obc2.xml')
        # Set up the system
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=2*unit.nanometers, constraints=app.HBonds,
                                        hydrogenMass=1.5*unit.amu, implicitSolventKappa=1.0/unit.nanometer) # Use HMR
    
    # Configure the integrator to use for the simulation
    temperature = 298.15*unit.kelvin
    friction = 1.0/unit.picoseconds
    # Use 4 fs for HMR
    timestep = 4.0*unit.femtoseconds
    integrator = LangevinMiddleIntegrator(temperature, friction, timestep)
    integrator.setConstraintTolerance(0.00001)  # Set the constraint error tolerance

    # Initialize the simulation
    simulation = app.Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('CUDA'))

    # Set the initial positions
    simulation.context.setPositions(pdb.positions)

    # Minimize the energy
    simulation.minimizeEnergy()

    # Identify non-water atoms for the HDF5Reporter
    non_water_atoms = [atom.index for atom in pdb.topology.atoms() if atom.residue.name not in ['HOH', 'WAT']]
    freq = 12500
    hdf5_reporter = md.reporters.HDF5Reporter((output_folder / 'output.h5').as_posix(), reportInterval=freq, atomSubset=non_water_atoms)
    statedata_reporter = app.StateDataReporter((output_folder / 'output.csv').as_posix(), reportInterval=freq, step=True, potentialEnergy=True, temperature=True,
                                               time = True, density = True, speed = True)
    checkput_reporter = app.CheckpointReporter((output_folder / 'checkput.chk').as_posix(), reportInterval=freq)
    simulation.reporters.append(hdf5_reporter)
    simulation.reporters.append(statedata_reporter)
    simulation.reporters.append(checkput_reporter)

    # Equilibration: running a short simulation to relax the structure
    equilibration_steps = np.round((simulation_time*unit.nanoseconds)/(4*unit.femtoseconds))  # Number of steps to equilibrate
    
    # Time the simulation for a short run of 1000 steps
    start_time = time.time()
    short_run = 12500*4
    simulation.step(short_run)
    end_time = time.time()

    # Calculate the elapsed time for the short run
    elapsed_time = end_time - start_time

    # Estimate the total wall time for the equilibration run based on the short run
    estimated_total_time = (elapsed_time / short_run) * (equilibration_steps-short_run)
    
    hours = estimated_total_time/60/60
    print(f"Takes {hours} hours.")
    sys.stdout.flush()
    
    simulation.step(equilibration_steps-short_run)

    # Don't forget to close the reporter when done to ensure data is properly saved
    hdf5_reporter.close()

if __name__ == '__main__':
    main()
