"""
cefe_lab.lattice -- spatial lattices and graph Laplacians.

Conventions
-----------
A lattice is a set of sites on a regular Cartesian grid inside a ball (3D)
or disk (2D).  `graph_laplacian` returns the POSITIVE graph Laplacian
( diag = sum of link weights, off-diag = -weight ), so the Klein-Gordon
operator is  W = L + m^2  (positive definite).

Link weights are 1/l^2 with l the link length; pass a `metric` weight
function w(i, j, coords) to curve the space (see
`schwarzschild_proper_distance`).
"""

import numpy as np


class Lattice:
    """A set of lattice sites with coordinates and grid indices."""
    def __init__(self, coords, indices, spacing, dim):
        self.coords = np.asarray(coords, dtype=float)   # (N, dim)
        self.indices = indices                           # list of grid tuples
        self.spacing = float(spacing)
        self.dim = dim

    def __len__(self):
        return len(self.coords)

    @property
    def radii(self):
        return np.linalg.norm(self.coords, axis=1)


def ball_lattice(radius, spacing):
    """3D Cartesian lattice inside the ball r <= radius."""
    ns = int(np.ceil(radius / spacing))
    coords, indices = [], []
    for ix in range(-ns, ns + 1):
        for iy in range(-ns, ns + 1):
            for iz in range(-ns, ns + 1):
                x, y, z = ix * spacing, iy * spacing, iz * spacing
                if np.sqrt(x * x + y * y + z * z) <= radius:
                    coords.append((x, y, z))
                    indices.append((ix, iy, iz))
    return Lattice(coords, indices, spacing, 3)


def disk_lattice(radius, spacing):
    """2D Cartesian lattice inside the disk r <= radius."""
    ns = int(np.ceil(radius / spacing))
    coords, indices = [], []
    for ix in range(-ns, ns + 1):
        for iy in range(-ns, ns + 1):
            x, y = ix * spacing, iy * spacing
            if np.sqrt(x * x + y * y) <= radius:
                coords.append((x, y, 0.0))
                indices.append((ix, iy))
    return Lattice(coords, indices, spacing, 2)


def schwarzschild_proper_distance(r1, r2, Rs, n=8):
    """Proper length of a radial link, int dr / sqrt(1 - Rs/r), Simpson."""
    xs = np.linspace(r1, r2, n + 1)
    fx = 1.0 / np.sqrt(1.0 - Rs / xs)
    h = (r2 - r1) / n
    return h / 3.0 * (fx[0] + fx[-1] + 4 * fx[1:-1:2].sum() + 2 * fx[2:-1:2].sum())


def graph_laplacian(lat, weight=None):
    """Positive graph Laplacian of the 6/4-neighbor Cartesian lattice.

    weight(i, j, coords) -> optional link-weight override; default 1/a^2.
    A Schwarzschild weight for radial links, e.g.:
        w = a / proper_distance  *  (1/a**2)
    """
    n = len(lat)
    a = lat.spacing
    index = {idx: i for i, idx in enumerate(lat.indices)}
    L = np.zeros((n, n))
    base = 1.0 / a ** 2
    steps = [(1, 0), (0, 1)] if lat.dim == 2 else [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    radii = lat.radii
    for i, idx in enumerate(lat.indices):
        for d in steps:
            jdx = tuple(idx[k] + d[k] for k in range(len(d)))
            j = index.get(jdx)
            if j is None or j < i:
                continue
            w = base if weight is None else weight(i, j, lat.coords, radii)
            L[i, j] -= w
            L[j, i] -= w
            L[i, i] += w
            L[j, j] += w
    return L
