import numpy as np

def __lawofCosines(a, b, c):
    """ Computes the angle between a and b in triangle ABC using the law of cosines

    Parameters:
        a (float): the length of side a
        b (float): the length of side b
        c (float): the length of side c

    Returns:
        The angle between a and b as a float
    """
    
    return np.arccos((a**2 + b**2 - c**2)/(2*a*b))

def __rodriguesRotation(v, theta, k):
    """ Rotates vector v by theta about unit vector k

    Parameters:
        v (numpy array): Vector to be rotated
        theta (float): Angle by which to rotate v
        k (numpy array): Unit vector about which to rotate
        
    Returns:
        Rotated v as a numpy array
    """
    
    return v*np.cos(theta) + np.cross(k, v)*np.sin(theta) + k*np.dot(k, v)*(1 - np.cos(theta))


def azct(origin, alpha1, alpha2, la, lb, lc, gamma, ssa=True, lb0=None):
    """ Uses inputs to generate a pair of asymmetric zipper-coupled tube segments
        with their smooth sheet attachments

    Parameters:
        origin (numpy array): Vector containing the origin of the segment
        alpha1 (float): First design angle in base unit: a degree-four vertex
        alpha1 (float): First design angle in base unit: a degree-four vertex
        la (float): Length of crease vector a
        lb (float): Length of crease vector b
        lc (float): Length of crease vector c
        gamma (float): Driving angle which determines how folded the structure is
        lb0 (float): Length of crease vector b in Z0 if using non-uniform extensions

    Returns:
        A (60, 6) numpy array containing all the edges in the component
    """
    

    """ ASYMMETRIC ZIPPER-COUPLED TUBES SETUP """
    # Derive dependent angles and length from inputs
    alpha3 = np.pi - alpha1
    alpha4 = np.pi - alpha2
    ld = lb*(np.sin(alpha1)/np.sin(alpha4))

    # Define constants for deriving important values
    ax = (np.sin(alpha1)*np.cos(alpha4) - np.sin(alpha4)*np.cos(alpha1))
    bx = (np.sin(alpha1)*np.cos(alpha3) - np.sin(alpha4)*np.cos(alpha2))
    ay = (np.sin(alpha1)*np.cos(alpha4) + np.sin(alpha4)*np.cos(alpha1))
    by = (np.sin(alpha1)*np.cos(alpha3) + np.sin(alpha4)*np.cos(alpha2))
    k = 2*np.sin(alpha1)*np.sin(alpha4)

    # Define the four crease vectors
    def a(g=gamma): return la*np.array([0, 0, -1])

    def b(g=gamma): return lb*np.array([-np.sin(alpha1)*np.cos(g),
                                   np.sin(alpha1)*np.sin(g),
                                   -np.cos(alpha1)])
    def c(g=gamma):
        
        # Use the above constants to compute the components of c
        c3 = -(ax*bx*(np.sin(g)**2) + ay*by*(np.cos(g)**2) +\
              (k**2)*(np.sin(g)**2)*(np.cos(g)**2))/\
              ((ax**2)*(np.sin(g)**2) + (ay**2)*(np.cos(g)**2) +\
               (k**2)*(np.sin(g)**2)*(np.cos(g)**2))
        c2 = (ay*c3 + by)/(k*np.sin(g))
        c1 = (ax*c3 + bx)/(k*np.cos(g))

        return lc*np.array([c1, c2, c3])
    
    def d(g=gamma): return ld*np.array([np.sin(alpha4)*np.cos(g),
                                   np.sin(alpha4)*np.sin(g),
                                   -np.cos(alpha4)])

    # Calculate the amount by which the second tube must be shifted
    if lb0 is not None:
        ld0 = lb0*(np.sin(alpha1)/np.sin(alpha4))
        ls = ((2*lb - lb0)*np.cos(alpha1) + (2*ld - ld0)*np.cos(alpha4))/2
    else:
        ls = (lb*np.cos(alpha1) + ld*np.cos(alpha4))/2

    # Define the shift vector s
    def s(g=gamma): return (1 + ls)*a(g) + c(g)

    # Define rotated crease vectors
    def abar(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@a(g)
    def bbar(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@b(g)
    def cbar(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@c(g)
    def dbar(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@d(g)

    # Define transformation for rotating the structure to be parallel to the xy-plane
    def phi(g=gamma): return np.arctan(-ax/(k*np.cos(g)))
    def Phi(g=gamma): return np.array([[np.cos(phi(g)), 0, np.sin(phi(g))],
                                       [0, 1, 0],
                                       [-np.sin(phi(g)), 0, np.cos(phi(g))]])

    """ SMOOTH SHEET SETUP """
    if ssa:
        
        # Derive gamma0--the value of gamma at which the tube reaches its ideal state
        gamma0 = np.arccos(np.sqrt((ay*by - ax**2 + k**2 + \
                                 np.sqrt((ay*by + ax**2 + k**2)**2 - \
                                         4*(ay**2)*ax*bx))/(2*(k**2))))

        # Define qbar and qperp
        def qbar(g): return dbar(g) - bbar(g)
        def qperp(g): return qbar(g) - np.dot(qbar(g), cbar(g))*cbar(g)/(lc**2)

        # Define Delta and delta:
        # important distances at gamma = gamma0 and gamma = pi/2
        Delta = np.linalg.norm(qperp(gamma0))
        delta = (ld*np.sin(alpha1) + lb*np.sin(alpha4) - Delta)/2

        # Define lengths for u and v
        lu1 = (ld*np.sin(alpha1) - lb*np.sin(alpha4) + Delta)/2
        lv1 = (lb*np.sin(alpha4) - ld*np.sin(alpha1) + Delta)/2
        lu2 = lu1/np.tan(alpha1)

        # Define the angle of rotation for u1
        def lam(g): return __lawofCosines(np.linalg.norm(qperp(g)), lu1, lv1)

        # Define the components of ubar
        def u1(g): return __rodriguesRotation(qperp(g)/np.linalg.norm(qperp(g)),
                                              lam(g), -cbar(g)/lc)
        def u2(g): return -cbar(g)/lc

        # Define ubar and vbar
        def ubar(g=gamma): return lu1*u1(g) + lu2*u2(g)
        def vbar(g=gamma): return ubar(g) - qbar(g)

        # Derive angles for computing candidate length lw1
        beta1 = np.arccos(np.dot(-cbar(gamma0), vbar(gamma0))/\
                          (lc*np.linalg.norm(vbar(gamma0))))
        beta2 = np.arccos(np.dot(-cbar(gamma0), d(gamma0) - ls*a(gamma0)/la)/\
                          (lc*np.linalg.norm(d(gamma0) - ls*a(gamma0)/la)))

        # Get the length of w
        lw1 = lc + lv1/np.tan(beta1) - lv1/np.tan(beta2)
        lw2 = lc - delta/np.tan(alpha1) - delta/np.tan(alpha2 - alpha1)
        lw = min(lw1, lw2)

        # Define wbar
        def wbar(g=gamma): return -lw*cbar(g)/lc

        # Define rotated crease vectors for smooth sheet attachment on first tube
        def u(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@ubar(g)
        def w(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@wbar(g)


    """ DEFINE VERTICES AND EDGES """

    # Tube 1
    E0 = d()
    E1 = a() + d()
    E2 = a() + c() + d()
    E3 = c() + d()
    
    F0 = np.array([0, 0, 0])
    F1 = a()
    F2 = a() + c()
    F3 = c()

    G0 = b()
    G1 = a() + b()
    G2 = a() + b() + c()
    G3 = b() + c()

    # Tube 2
    E0p = dbar() + s()
    E1p = abar() + dbar() + s()
    E2p = abar() + cbar() + dbar() + s()
    E3p = cbar() + dbar() + s()
    
    F0p = s()
    F1p = abar() + s()
    F2p = abar() + cbar() + s()
    F3p = cbar() + s()

    G0p = bbar() + s()
    G1p = abar() + bbar() + s()
    G2p = abar() + bbar() + cbar() + s()
    G3p = bbar() + cbar() + s()

    if ssa:
        # Smooth Sheet Attachment 1
        H0 = u()
        H1 = a() + c() + d() - u() + w()
        H2 = a() + c() + d() - u()
        H3 = u() - w()

        # Smooth Sheet Attachment 2
        H0p = s() + ubar()
        H1p = abar() + cbar() + dbar() + s() - ubar() + wbar()
        H2p = abar() + cbar() + dbar() + s() - ubar()
        H3p = s() + ubar() - wbar()

        # Points on adjacent tubes
        F0_1 = d() - b()
        F3_1 = c() + d() - b()
        F0p_n1 = s() + dbar() - bbar()
        F3p_n1 = s() + cbar() + dbar() - bbar()
    
    # Define all the edges in the structure
    edges = [(E0, E1), (E1, E2), (E2, E3), (E3, E0), # Tube 1
             (F0, F1), (F1, F2), (F2, F3), (F3, F0),
             (G0, G1), (G1, G2), (G2, G3), (G3, G0),
             (E0, F0), (E1, F1), (E2, F2), (E3, F3),
             (F0, G0), (F1, G1), (F2, G2), (F3, G3),
             
             (E0p, E1p), (E1p, E2p), (E2p, E3p), (E3p, E0p), # Tube 2
             (F0p, F1p), (F1p, F2p), (F2p, F3p), (F3p, F0p),
             (G0p, G1p), (G1p, G2p), (G2p, G3p), (G3p, G0p),
             (E0p, F0p), (E1p, F1p), (E2p, F2p), (E3p, F3p),
             (F0p, G0p), (F1p, G1p), (F2p, G2p), (F3p, G3p)]

    if ssa:
        sheet_edges = [(E1, H1), (E2, H2), (G1, H1), (G2, H2), # Sheet 1
                       (F0, H0), (F3, H3), (F0_1, H0), (F3_1, H3),
                       (H0, H3), (H1, H2),
                       
                       (E1p, H1p), (E2p, H2p), (G1p, H1p), (G2p, H2p), # Sheet 2
                       (F0p, H0p), (F3p, H3p), (F0p_n1, H0p), (F3p_n1, H3p),
                       (H0p, H3p), (H1p, H2p)]
        edges = np.concatenate((edges, sheet_edges))

    # Compose all the edges into a single list for plotting
    # Apply the transformation Phi
    component = np.empty((0, 6))
    for edge in edges:
        component = np.vstack((component, np.concatenate((Phi()@edge[0] + Phi()@origin,
                                                          Phi()@edge[1] - Phi()@edge[0]))))

    return component
