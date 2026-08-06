#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh

# Ensure cefe_py is importable
try:
    import cefe_py as ce
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
    import cefe_py as ce

def build_python_hamiltonian(points_3d, spacing, alpha=10.0):
    n = len(points_3d)
    print(f"Building Hamiltonian for {n} spatial points...")
    
    # Extract coordinates
    coords = np.array([[p.x, p.y, p.z] for p in points_3d])
    
    # KDTree for fast neighbor lookup
    tree = cKDTree(coords)
    
    H = lil_matrix((n, n), dtype=float)
    
    # Find neighbors within distance slightly larger than spacing
    search_radius = spacing * 1.1
    neighbors = tree.query_ball_point(coords, search_radius)
    
    # Laplacian coefficients
    inv_a2 = 1.0 / (spacing * spacing)
    
    for i in range(n):
        r = np.linalg.norm(coords[i])
        
        # Entropic gravity inverse-distance potential: V(r) = - alpha / r
        # Add a small epsilon to avoid singularity at origin
        v_r = -alpha / (r + 0.1 * spacing)
        
        diag_val = v_r
        num_neighbors = 0
        
        for j in neighbors[i]:
            if i != j:
                # -0.5 * Laplacian for kinetic energy (hbar=1, m=1)
                # L_ij = inv_a2, H_ij = -0.5 * inv_a2
                H[i, j] = -0.5 * inv_a2
                num_neighbors += 1
                
        # Diagonal term of Laplacian is -num_neighbors * inv_a2
        # So H_ii = -0.5 * (-num_neighbors * inv_a2) + V
        H[i, i] = 0.5 * num_neighbors * inv_a2 + v_r
        
    return H.tocsr(), coords

def run_entropic_atom():
    print("Initializing Causal Diamond for Entropic Atom...")
    grid = ce.CausalDiamondGrid(radius=10.0, spacing=0.5, metric=ce.MinkowskiMetric())
    
    # Get all points
    all_points = grid.get_points()
    # Filter points on t=0 slice
    points_3d = [p for p in all_points if abs(p.t) < 1e-9]
    
    # Build Hamiltonian with an emergent entropic potential alpha=5.0
    H, coords = build_python_hamiltonian(points_3d, spacing=0.5, alpha=5.0)
    
    print("Solving for emergent orbital shapes (lowest eigenstates)...")
    # Find 4 smallest algebraic eigenvalues
    evals, evecs = eigsh(H, k=4, which='SA')
    
    print(f"Ground state energy: {evals[0]:.4f}")
    print(f"Excited state energies: {evals[1]:.4f}, {evals[2]:.4f}, {evals[3]:.4f}")
    
    # Visualize
    # We want to plot the 2D cross section at z ~ 0
    z_mask = np.abs(coords[:, 2]) < 0.25
    idx_2d = np.where(z_mask)[0]
    x_2d = coords[idx_2d, 0]
    y_2d = coords[idx_2d, 1]
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 10), facecolor='#0a0a0f')
    fig.suptitle('Emergent Entropic Atom Orbitals on Causal Grid', color='white', fontsize=16)
    
    for i, ax in enumerate(axes.flatten()):
        psi = evecs[:, i]
        prob = np.abs(psi[idx_2d])**2
        
        # Normalize probability for better coloring
        if np.max(prob) > 0:
            prob = prob / np.max(prob)
            
        ax.set_facecolor('#111118')
        sc = ax.scatter(x_2d, y_2d, c=prob, cmap='magma', s=25, alpha=0.9, vmin=0, vmax=1)
        ax.set_title(f'State n={i} (E={evals[i]:.2f})', color='#ff9900')
        ax.set_aspect('equal')
        ax.set_xlim(-10, 10)
        ax.set_ylim(-10, 10)
        for spine in ax.spines.values():
            spine.set_color('#333333')
        ax.tick_params(colors='#666666')
        
    plt.tight_layout()
    plt.savefig('entropic_atom_orbitals.png', facecolor='#0a0a0f', dpi=150)
    print("Saved emergent shapes to entropic_atom_orbitals.png")

if __name__ == "__main__":
    run_entropic_atom()
