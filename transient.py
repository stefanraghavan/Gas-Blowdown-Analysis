import matplotlib.pyplot as plt
import numpy as np

# Constants and parameters
Diameter = 0.0417  # Pipe diameter in feet
Length = 10.0  # Pipe length in feet
GridPoints = 101  # Number of grid points along the pipe
TimeSteps = 1000  # Number of time steps
DeltaX = Length / (GridPoints - 1)
DeltaT = 0.01  # Time step size (adjust as needed)
FluidDensity = 998.2  # Water density in kg/m^3
Viscosity = 0.001  # Water dynamic viscosity in kg/(m·s)

# Initialize arrays for pressure, velocity, and density
Pressure = np.zeros(GridPoints)
Velocity = np.zeros(GridPoints)
Density = np.full(GridPoints, FluidDensity)
initial_pressure = 100.0  # Example initial pressure
Pressure.fill(initial_pressure)

# Define initial conditions for velocity
max_velocity = 1.0  # Maximum velocity at the pipe center (you can adjust this value)
for i in range(GridPoints):
    # Parabolic velocity profile from zero at the walls to max_velocity at the center
    Velocity[i] = max_velocity * (1.0 - (2.0 * (i / (GridPoints - 1) - 0.5))**2)

# Time-stepping loop
for t in range(TimeSteps):
    # Update pressure and velocity fields using FDM discretization
    for i in range(1, GridPoints - 1):
        # Implement the discretized equations here
        Velocity[i] = Velocity[i] - (DeltaT / DeltaX) * (Density[i] - Density[i-1])
        Pressure[i] = Pressure[i] - (DeltaT / DeltaX) * (Density[i] * Velocity[i+1] - Density[i] * Velocity[i]) + (Viscosity * DeltaT / DeltaX**2) * (Velocity[i+1] - 2 * Velocity[i] + Velocity[i-1])
    
    # Apply boundary conditions (you should implement these)
    # For example, set the inlet and outlet boundary conditions here
    # Apply outlet boundary condition (fixed pressure outlet)
    outlet_pressure = 90.0  # Specify the desired outlet pressure (adjust as needed)
    Pressure[-1] = outlet_pressure
    
    # Track pressure at the end of the pipe
    pressure_at_end = Pressure[-1]

    # Print or save results at each time step
    if t % 100 == 0:
        print(f"Time step {t}: Pressure at the end = {pressure_at_end}")

# Plot pressure distribution at the end of the simulation
plt.plot(Pressure)
plt.xlabel("Grid Points")
plt.ylabel("Pressure")
plt.title("Pressure Distribution along the Pipe")
plt.show()
