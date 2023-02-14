import numpy as np

# Constants
eps = 1e-12

# Node class for the tree


class Node:
    def __init__(self, x, y, level, parent=None):
        self.x = x
        self.y = y
        self.children = []
        self.level = level
        self.parent = parent
        self.is_leaf = True
        self.charge = 0.0
        self.mx = 0.0
        self.my = 0.0

    def add_child(self, child):
        self.children.append(child)
        self.is_leaf = False

# FMM class


class FMM:
    def __init__(self, charges, x, y):
        self.charges = charges
        self.x = x
        self.y = y
        self.tree = None
        self.max_level = None
        self.build_tree()
        self.compute_multipole_moments()
        self.compute_forces()

    def build_tree(self):
        # Find the bounding box
        xmin = self.x.min()
        xmax = self.x.max()
        ymin = self.y.min()
        ymax = self.y.max()

        # Determine the maximum level
        self.max_level = int(np.ceil(np.log2(max(xmax-xmin, ymax-ymin)/eps)))

        # Create the root node
        self.tree = Node((xmin+xmax)/2, (ymin+ymax)/2, 0)

        # Add the charges to the tree
        for i in range(len(self.charges)):
            node = self.tree
            for level in range(self.max_level):
                # Determine which quadrant the charge belongs to
                if self.x[i] < node.x:
                    if self.y[i] < node.y:
                        quadrant = 0
                    else:
                        quadrant = 1
                else:
                    if self.y[i] < node.y:
                        quadrant = 2
                    else:
                        quadrant = 3

                # If there is no child in the quadrant, create one
                if len(node.children) <= quadrant:
                    xmin = node.x if quadrant in [0, 1] else node.x
                    xmax = node.x if quadrant in [2, 3] else node.x
                    ymin = node.y if quadrant in [0, 2] else node.y
                    ymax = node.y if quadrant in [1, 3] else node.y
                    node.add_child(
                        Node((xmin+xmax)/2, (ymin+ymax)/2, level+1, node))

                # Move to the child node
                node = node.children[quadrant]

            # Add the charge to the leaf node
            node.charge += self.charges[i]

    def compute_multipole_moments(self):
        for node in self.
