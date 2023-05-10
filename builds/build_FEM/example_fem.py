import numpy as np

# Define the properties of the truss system
E = 200000  # Young's modulus in MPa
A = 0.01  # Cross-sectional area in m^2

# Define the nodes and elements
nodes = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])  # Nodes coordinates
# Element connectivity
elements = np.array([[0, 1], [0, 2], [1, 3], [2, 3], [0, 3]])

# Define the loads and supports
loads = np.array([[0, 0], [0, 0], [0, -1], [0, 0]])  # Loads on the nodes
# Supports at the nodes (1 if supported, 0 if free)
supports = np.array([[1, 1], [1, 0], [0, 0], [0, 0]])

# Initialize global stiffness matrix and force vector
K_global = np.zeros((2*len(nodes), 2*len(nodes)))
F_global = loads.flatten()

# Loop over elements to assemble global stiffness matrix
for element in elements:
    node1, node2 = nodes[element[0]], nodes[element[1]]
    L = np.linalg.norm(node2 - node1)
    direction_cosine = (node2 - node1) / L
    transformation_matrix = np.outer(direction_cosine, direction_cosine)

    # Local stiffness matrix for 2D truss element
    K_local = (E * A / L) * np.array([
        [transformation_matrix[0, 0], transformation_matrix[0, 1], -
            transformation_matrix[0, 0], -transformation_matrix[0, 1]],
        [transformation_matrix[1, 0], transformation_matrix[1, 1], -
            transformation_matrix[1, 0], -transformation_matrix[1, 1]],
        [-transformation_matrix[0, 0], -transformation_matrix[0, 1],
            transformation_matrix[0, 0], transformation_matrix[0, 1]],
        [-transformation_matrix[1, 0], -transformation_matrix[1, 1],
            transformation_matrix[1, 0], transformation_matrix[1, 1]]
    ])

    assembly_indices = np.array(
        [2*element[0], 2*element[0]+1, 2*element[1], 2*element[1]+1])
    K_global[np.ix_(assembly_indices, assembly_indices)] += K_local

# Apply boundary conditions
bc_indices = np.where(supports.flatten() == 1)[0]
K_global = np.delete(K_global, bc_indices, axis=0)
K_global = np.delete(K_global, bc_indices, axis=1)
F_global = np.delete(F_global, bc_indices)

# Solve for displacements
displacements = np.linalg.solve(K_global, F_global)

print(displacements)
