# src/solver.py
import numpy as np

def compute_flux(u_state, direction='x', gamma=1.4):
    """Compute Euler fluxes in x or y direction."""
    rho = u_state[0]
    u_vel = u_state[1] / rho
    v_vel = u_state[2] / rho
    E = u_state[3]
    p = (gamma - 1) * (E - 0.5 * rho * (u_vel**2 + v_vel**2))
    
    flux = np.zeros_like(u_state)
    if direction == 'x':
        flux[0] = rho * u_vel
        flux[1] = rho * u_vel**2 + p
        flux[2] = rho * u_vel * v_vel
        flux[3] = u_vel * (E + p)
    else:  # y direction
        flux[0] = rho * v_vel
        flux[1] = rho * u_vel * v_vel
        flux[2] = rho * v_vel**2 + p
        flux[3] = v_vel * (E + p)
    return flux

def hllc_flux(u_L, u_R, direction='x', gamma=1.4):
    """Compute HLLC Riemann solver flux for a single interface."""
    if u_L.shape != (4,) or u_R.shape != (4,):
        raise ValueError(f"Invalid input shape: u_L={u_L.shape}, u_R={u_R.shape}, expected (4,)")
    
    # Left state
    rho_L = u_L[0]
    u_vel_L = u_L[1] / rho_L
    v_vel_L = u_L[2] / rho_L
    E_L = u_L[3]
    p_L = (gamma - 1) * (E_L - 0.5 * rho_L * (u_vel_L**2 + v_vel_L**2))
    a_L = np.sqrt(max(gamma * p_L / rho_L, 0.0))
    
    # Right state
    rho_R = u_R[0]
    u_vel_R = u_R[1] / rho_R
    v_vel_R = u_R[2] / rho_R
    E_R = u_R[3]
    p_R = (gamma - 1) * (E_R - 0.5 * rho_R * (u_vel_R**2 + v_vel_R**2))
    a_R = np.sqrt(max(gamma * p_R / rho_R, 0.0))
    
    # Normal/tangential velocities
    if direction == 'x':
        u_n_L, u_t_L = u_vel_L, v_vel_L
        u_n_R, u_t_R = u_vel_R, v_vel_R
    else:
        u_n_L, u_t_L = v_vel_L, u_vel_L
        u_n_R, u_t_R = v_vel_R, u_vel_R
    
    # Estimate p_star and u_star (Toro-style)
    denom = rho_L * a_L + rho_R * a_R
    if denom == 0:
        fL = compute_flux(u_L, direction, gamma)
        fR = compute_flux(u_R, direction, gamma)
        return 0.5 * (fL + fR) - 0.5 * max(a_L, a_R) * (u_R - u_L)
    
    p_star = (rho_L * a_L * p_R + rho_R * a_R * p_L + rho_L * a_L * rho_R * a_R * (u_n_L - u_n_R)) / denom
    u_star = (u_n_L * rho_L * a_L + u_n_R * rho_R * a_R + p_L - p_R) / denom
    
    # Wave speed estimates
    S_L = min(u_n_L - a_L, u_star)
    S_R = max(u_n_R + a_R, u_star)
    S_star = u_star
    
    # Flux region
    if S_L >= 0:
        return compute_flux(u_L, direction, gamma)
    if S_R <= 0:
        return compute_flux(u_R, direction, gamma)
    
    # Star states
    u_star_L = np.zeros(4)
    u_star_R = np.zeros(4)
    rho_star_L = rho_L * (S_L - u_n_L) / (S_L - S_star)
    rho_star_R = rho_R * (S_R - u_n_R) / (S_R - S_star)
    p_star_calc = p_L + rho_L * (S_L - u_n_L) * (S_star - u_n_L)
    
    if direction == 'x':
        u_star_L[0] = rho_star_L
        u_star_L[1] = rho_star_L * S_star
        u_star_L[2] = rho_star_L * u_t_L
        u_star_L[3] = p_star_calc / (gamma - 1) + 0.5 * rho_star_L * (S_star**2 + u_t_L**2)
        u_star_R[0] = rho_star_R
        u_star_R[1] = rho_star_R * S_star
        u_star_R[2] = rho_star_R * u_t_R
        u_star_R[3] = p_star_calc / (gamma - 1) + 0.5 * rho_star_R * (S_star**2 + u_t_R**2)
    else:
        u_star_L[0] = rho_star_L
        u_star_L[1] = rho_star_L * u_t_L
        u_star_L[2] = rho_star_L * S_star
        u_star_L[3] = p_star_calc / (gamma - 1) + 0.5 * rho_star_L * (u_t_L**2 + S_star**2)
        u_star_R[0] = rho_star_R
        u_star_R[1] = rho_star_R * u_t_R
        u_star_R[2] = rho_star_R * S_star
        u_star_R[3] = p_star_calc / (gamma - 1) + 0.5 * rho_star_R * (u_t_R**2 + S_star**2)
    
    if S_star >= 0:
        return compute_flux(u_L, direction, gamma) + S_L * (u_star_L - u_L)
    return compute_flux(u_R, direction, gamma) + S_R * (u_star_R - u_R)