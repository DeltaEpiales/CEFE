import numpy as np
from cefe_py import cefe_core as ce

def test_dynamic_tensor_backreaction():
    print("--- Testing Dynamic Tensor Engine (True Backreaction) ---", flush=True)
    
    # 1. Initialize CausalDiamondGrid (using Minkowski as base)
    radius = 3.0
    spacing = 0.5
    grid = ce.CausalDiamondGrid(radius, spacing, ce.MinkowskiMetric())
    
    # 2. Initialize DynamicTensorEngine
    print(ce.__file__, flush=True)
    print(dir(ce), flush=True)
    dte = ce.DynamicTensorEngine(grid)
    
    # 3. Create a mock scalar field configuration (dense at origin)
    nodes = grid.get_points()
    N = len(nodes)
    field_phi = np.zeros(N)
    for i, node in enumerate(nodes):
        # A gaussian wavepacket centered at origin
        r = np.sqrt(node.x**2 + node.y**2 + node.z**2)
        field_phi[i] = np.exp(-r**2)
        
    print(f"Grid initialized with {N} topological nodes.", flush=True)
    
    # Get initial r-coordinate of a node near the origin (not exactly 0)
    test_node_idx = None
    for i, node in enumerate(nodes):
        node_r = np.sqrt(node.x**2 + node.y**2 + node.z**2)
        if node_r > 0.1 and node_r < 1.0:
            test_node_idx = i
            break
            
    if test_node_idx is not None:
        initial_node = nodes[test_node_idx]
        initial_r = np.sqrt(initial_node.x**2 + initial_node.y**2 + initial_node.z**2)
        print(f"Initial r-coordinate of node {test_node_idx}: {initial_r:.6f}", flush=True)
        
        # 4. Compute Stress-Energy Tensor and warp geometry
        mass = 1.0
        time_step = 0.1
        G_constant = 1.0
        
        print("Computing Stress-Energy Tensor...", flush=True)
        dte.compute_stress_energy_tensor(field_phi, mass)
        print("Applying Linearized Einstein Update...", flush=True)
        dte.linearized_einstein_update(time_step, G_constant)
        print("Finished Einstein Update.", flush=True)
        
        # 5. Verify the geometry warped
        warped_nodes = grid.get_points()
        final_node = warped_nodes[test_node_idx]
        final_r = np.sqrt(final_node.x**2 + final_node.y**2 + final_node.z**2)
        print(f"Final r-coordinate of node {test_node_idx}: {final_r:.6f}", flush=True)
        
        if final_r < initial_r:
            print("SUCCESS: Geometry correctly warped inwards due to Energy Density (Gravity is attractive).", flush=True)
        else:
            print("FAILURE: Geometry did not warp as expected.", flush=True)
    else:
        print("Could not find suitable test node.", flush=True)

if __name__ == "__main__":
    test_dynamic_tensor_backreaction()
