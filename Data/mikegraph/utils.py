"""
Utils module for mikegraph
"""
import numpy as np

def calculate_full_flow(diameter, slope, material, resolution=1e-6):
    """
    Calculate full-pipe discharge capacity using the Colebrook-White equation.

    Determines the maximum flow capacity of a circular pipe running full using
    iterative solution of the Colebrook-White friction factor equation.

    Parameters:
        diameter (float): Internal pipe diameter in meters
        slope (float): Hydraulic gradient (head loss per unit length), dimensionless
        material (str): Pipe material type. Materials starting with 'p' (e.g., 'plastic', 
            'PVC') use roughness k=0.001m. All others use k=0.0015m.
        resolution (float, optional): Convergence tolerance for iterative solution.
            Defaults to 1e-6.

    Returns:
        float or None: Full-pipe discharge in m³/s. Returns None if solution 
            doesn't converge within the specified resolution.

    Notes:
        - Uses kinematic viscosity of 1.3e-6 m²/s (water at ~10°C)
        - Employs bisection method to find velocity that satisfies flow equations
        - Iteratively solves Colebrook-White equation for friction factor
        - Hydraulic radius R = D/4 for circular pipes running full
        - Flow area = π(D/2)² for circular cross-section
    Examples:
        >>> Q = calculate_full_flow(0.30, 0.005, "PVC")
    """
    if material[0].lower() == "p":
        k = 0.001
    else:
        k = 0.0015

    g = 9.82  # m2/s
    kinematicViscosity = 0.0000013  # m2/s
    R = diameter / 4.0

    Vmin = 0.001
    Vmax = 500
    I = 1E+40
    while abs(slope - I) > resolution and Vmax - Vmin > resolution:
        V = 10 ** ((np.log10(Vmax) + np.log10(Vmin)) / 2.0)

        Re = V * R / kinematicViscosity
        f = 0.01  # Initial Guess
        fOld = 1000
        while abs(f - fOld) > resolution:
            fOld = f
            f = 2 / (6.4 - 2.45 * np.log(k / R + 4.7 / (Re * np.sqrt(f)))) ** 2
        I = f * (V ** 2 / (2 * g * R))
        if slope > I:
            Vmin = V
        else:
            Vmax = V
    if Vmax - Vmin <= resolution:
        return None
    else:
        return V * (diameter / 2.0) ** 2 * np.pi
