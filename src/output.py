# src/output.py
import numpy as np
import matplotlib.pyplot as plt
import os

def save_output(u, t, X, Y, gamma=1.4, results_dir="results",
                rho_range=None, uvel_range=None, vvel_range=None, p_range=None):
    """Save simulation data and generate contour plots with fixed color scales."""
    os.makedirs(results_dir, exist_ok=True)
    
    rho = u[0, :, :]
    u_vel = u[1, :, :] / u[0, :, :]
    v_vel = u[2, :, :] / u[0, :, :]
    p = (gamma - 1) * (u[3, :, :] - 0.5 * u[0, :, :] * (u_vel**2 + v_vel**2))
    
    # Save data
    np.savez(f"{results_dir}/fields_t{t:.3f}.npz", rho=rho, u_vel=u_vel, v_vel=v_vel, p=p)
    
    # Plot with fixed ranges
    plt.figure(figsize=(12, 10))
    
    plt.subplot(2, 2, 1)
    plt.contourf(X, Y, u_vel.T, cmap='inferno',
                 vmin=uvel_range[0] if uvel_range else None,
                 vmax=uvel_range[1] if uvel_range else None)
    plt.colorbar(label='x-Velocity')
    plt.title('x-Velocity')

    plt.subplot(2, 2, 2)
    plt.contourf(X, Y, rho.T, cmap='viridis',
                 vmin=rho_range[0] if rho_range else None,
                 vmax=rho_range[1] if rho_range else None)
    plt.colorbar(label='Density')
    plt.title('Density')
    
    plt.subplot(2, 2, 3)
    plt.contourf(X, Y, v_vel.T, cmap='plasma',
                 vmin=vvel_range[0] if vvel_range else None,
                 vmax=vvel_range[1] if vvel_range else None)
    plt.colorbar(label='y-Velocity')
    plt.title('y-Velocity')
    
    plt.subplot(2, 2, 4)
    plt.contourf(X, Y, p.T, cmap='magma',
                 vmin=p_range[0] if p_range else None,
                 vmax=p_range[1] if p_range else None)
    plt.colorbar(label='Pressure')
    plt.title('Pressure')
    
    plt.suptitle(f'2D Shock Tube at t = {t:.3f}')
    plt.tight_layout()
    plt.savefig(f"{results_dir}/figure_t{t:.3f}.png")
    plt.close()