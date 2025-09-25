# src/main.py
import numpy as np
from grid import setup_grid
from solver import compute_flux, hllc_flux
from output import save_output

# Simulation parameters
nx, ny = 200, 200  # Grid points
Lx, Ly = 1.0, 1.0  # Domain size
t_end = 0.2  # End time
CFL = 0.8  # CFL number
output_interval = 0.01  # Output frequency
gamma = 1.4  # Specific heat ratio

def main():
    # Setup grid and initial conditions
    x, y, X, Y, u, dx, dy, gamma, rho_range, uvel_range, vvel_range, p_range = setup_grid(nx, ny, Lx, Ly)
    
    # Time-stepping loop
    t = 0.0
    u_new = np.copy(u)
    next_output_time = 0.0
    
    while t < t_end:
        # Compute time step
        max_speed = 1e-16
        for i in range(nx):
            for j in range(ny):
                rho = u[0, i, j]
                u_vel = u[1, i, j] / rho
                v_vel = u[2, i, j] / rho
                p = (gamma - 1) * (u[3, i, j] - 0.5 * rho * (u_vel**2 + v_vel**2))
                a = np.sqrt(gamma * p / rho)
                max_speed = max(max_speed, abs(u_vel) + a, abs(v_vel) + a)
        
        dt = CFL * min(dx, dy) / max_speed
        if t + dt > t_end:
            dt = t_end - t
        
        # Compute fluxes
        flux_x = np.zeros((4, nx + 1, ny))
        flux_y = np.zeros((4, nx, ny + 1))
        
        for j in range(ny):
            for i in range(nx - 1):
                flux_x[:, i + 1, j] = hllc_flux(u[:, i, j], u[:, i + 1, j], 'x', gamma)
        
        for i in range(nx):
            for j in range(ny - 1):
                flux_y[:, i, j + 1] = hllc_flux(u[:, i, j], u[:, i, j + 1], 'y', gamma)
        
        # Boundary conditions (outflow)
        for j in range(ny):
            flux_x[:, 0, j] = compute_flux(u[:, 0, j], 'x', gamma)
            flux_x[:, -1, j] = compute_flux(u[:, -1, j], 'x', gamma)
        for i in range(nx):
            flux_y[:, i, 0] = compute_flux(u[:, i, 0], 'y', gamma)
            flux_y[:, i, -1] = compute_flux(u[:, i, -1], 'y', gamma)
        
        # Update conservative variables
        for i in range(nx):
            for j in range(ny):
                u_new[:, i, j] = u[:, i, j] - (dt / dx) * (flux_x[:, i + 1, j] - flux_x[:, i, j]) \
                                              - (dt / dy) * (flux_y[:, i, j + 1] - flux_y[:, i, j])
        u[:] = u_new
        t += dt
        
        # Output results
        if t >= next_output_time - 1e-10:
            save_output(u, t, X, Y, gamma)
            save_output(u, t, X, Y, gamma=gamma, results_dir="results",\
                            rho_range=rho_range, uvel_range=uvel_range, \
                                vvel_range=vvel_range, p_range=p_range)
            next_output_time += output_interval

if __name__ == "__main__":
    main()