# src/output.py
import numpy as np
import matplotlib.pyplot as plt
import os

def save_output(u, t, X, Y, gamma=1.4, results_dir="results"):
    """Save simulation data and generate contour plots."""
    os.makedirs(results_dir, exist_ok=True)
    
    rho = u[0, :, :]
    u_vel = u[1, :, :] / u[0, :, :]
    v_vel = u[2, :, :] / u[0, :, :]
    p = (gamma - 1) * (u[3, :, :] - 0.5 * u[0, :, :] * (u_vel**2 + v_vel**2))
    
    # Save data
    np.savez(f"{results_dir}/fields_t{t:.3f}.npz", rho=rho, u_vel=u_vel, v_vel=v_vel, p=p)
    
    # Generate plots
    plt.figure(figsize=(12, 10))
    plt.subplot(2, 2, 1)
    plt.contourf(X, Y, rho.T, cmap='viridis')
    plt.colorbar(label='Density')
    plt.title('Density')
    
    plt.subplot(2, 2, 2)
    plt.contourf(X, Y, u_vel.T, cmap='inferno')
    plt.colorbar(label='x-Velocity')
    plt.title('x-Velocity')
    
    plt.subplot(2, 2, 3)
    plt.contourf(X, Y, v_vel.T, cmap='plasma')
    plt.colorbar(label='y-Velocity')
    plt.title('y-Velocity')
    
    plt.subplot(2, 2, 4)
    plt.contourf(X, Y, p.T, cmap='magma')
    plt.colorbar(label='Pressure')
    plt.title('Pressure')
    
    plt.suptitle(f'2D Shock Tube at t = {t:.3f}')
    plt.tight_layout()
    plt.savefig(f"{results_dir}/figure_t{t:.3f}.png")
    plt.close()