# src/grid.py
import numpy as np

def setup_grid(nx=200, ny=200, Lx=1.0, Ly=1.0):
    """Set up 2D grid and initial conditions for shock tube."""
    dx, dy = Lx / nx, Ly / ny
    x = np.linspace(0.5 * dx, Lx - 0.5 * dx, nx)  # Cell centers
    y = np.linspace(0.5 * dy, Ly - 0.5 * dy, ny)
    X, Y = np.meshgrid(x, y, indexing='xy')
    
    rho_range = (min(rho_L, rho_R), max(rho_L, rho_R))
    uvel_range = (-0.5, 0.5)    # pick expected range
    vvel_range = (-0.5, 0.5)
    p_range = (min(p_L, p_R), max(p_L, p_R))

    # Conservative variables: [rho, rho*u, rho*v, E]
    u = np.zeros((4, nx, ny))
    
    # Initial conditions (left and right states)
    rho_L, u_L, v_L, p_L = 1.0, 0.0, 0.0, 1.0  # Left state
    rho_R, u_R, v_R, p_R = 0.125, 0.0, 0.0, 0.1  # Right state
    gamma = 1.4  # Specific heat ratio (air)
    
    for i in range(nx):
        for j in range(ny):
            if x[i] < 0.5:
                rho, uu, vv, p = rho_L, u_L, v_L, p_L  # Left side
            else:
                rho, uu, vv, p = rho_R, u_R, v_R, p_R  # Right side
            u[0, i, j] = rho
            u[1, i, j] = rho * uu
            u[2, i, j] = rho * vv
            u[3, i, j] = p / (gamma - 1) + 0.5 * rho * (uu**2 + vv**2)
    
    return x, y, X, Y, u, dx, dy, gamma, rho_range, uvel_range, vvel_range, p_range