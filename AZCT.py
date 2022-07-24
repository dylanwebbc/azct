import numpy as np
from scipy import linalg

def miuraori(e, f, g):
    """ Solves for the fourth vector in a degree-four vertex
        given the other three

    Parameters:
        a (numpy array): The first vector
        b (numpy array): The second vector
        c (numpy array): The third vector

    Returns:
        The fourth vector in the degree-four vertex as a numpy array
    """

    # Normalize vectors
    e = e/np.linalg.norm(e)
    f = f/np.linalg.norm(f)
    g = g/np.linalg.norm(g)

    # Get components
    e1, e2, e3 = e
    g1, g2, g3 = g

    # Get angles
    theta1 = np.pi - np.arccos(np.dot(f, g))
    theta2 = np.pi - np.arccos(np.dot(e, f))

    # Define linear components
    ax = ((e2/g2)*g3 - e3)/(e1 - g1*e2/g2)
    bx = (np.cos(theta1) - (e2/g2)*np.cos(theta2))/(e1 - g1*e2/g2)

    ay = ((g1/e1)*e3 - g3)/(g2 - e2*g1/e1)
    by = (np.cos(theta2) - (g1/e1)*np.cos(theta1))/(g2 - e2*g1/e1)

    # Define quadratic components
    a = ax**2 + ay**2 + 1
    b = 2*ax*bx + 2*ay*by
    c = bx**2 + by**2 - 1

    # Solve for d
    # Solve for d
    if np.allclose(b**2 - 4*a*c, 0):
        z = -b/(2*a)
    else:
        z = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    y = ay*z + by
    x = ax*z + bx

    return np.array([x, y, z])

def intersection(p1, p2, v1, v2):
    """ Solves for the length of v1 from p1 to the
        point of intersection with v2 starting at p2

    Parameters:
        p1 (numpy array): Point 1
        p2 (numpy array): Point 2
        v1 (numpy array): Vector 1
        v2 (numpy array): Vector 2

    Returns:
        The length of v1 as a float
    """

    A = np.array([[v1[0], -v2[0]], [v1[1], -v2[1]]])
    b = np.array([p2[0] - p1[0], p2[1] - p1[1]])

    return linalg.solve(A, b)[0]

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


def azct(origin, alpha1, alpha2, la, lb, lc, gamma, ssa=0, lb0=None):
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
        ls = ((2*lb - lb0)*np.cos(alpha1) + (2*ld - ld0)*np.cos(alpha4))/(2*la)
    else:
        ls = (lb*np.cos(alpha1) + ld*np.cos(alpha4))/(2*la)

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

        # Define Delta, delta and epsilon:
        # important distances at gamma = gamma0 and gamma = pi/2
        # and the distance to extend edges to remove gaps
        Delta = np.linalg.norm(qperp(gamma0))
        delta = (ld*np.sin(alpha1) + lb*np.sin(alpha4) - Delta)/2
        eps = (lb*np.sin(alpha1)*np.sin(gamma0)/2)/abs(c(gamma0)[1]/lc)

        # Define lengths for u and v
        lu1 = (ld*np.sin(alpha1) - lb*np.sin(alpha4) + Delta)/2
        lv1 = (lb*np.sin(alpha4) - ld*np.sin(alpha1) + Delta)/2
        if ssa == 1:
            lu2 = lu1/np.tan(alpha1)        
        elif ssa == 2:
            lu2 = (lu1/Delta)*np.linalg.norm(qbar(gamma0) - qperp(gamma0))
        lv2 = lv1*lu2/lu1

        # Define the angle of rotation for u1
        def lam(g): return __lawofCosines(np.linalg.norm(qperp(g)), lu1, lv1)

        # Define the components of ubar
        def u1(g): return __rodriguesRotation(qperp(g)/np.linalg.norm(qperp(g)),
                                              lam(g), -cbar(g)/lc)
        def u2(g): return -cbar(g)/lc

        # Define ubar, vbar and epsilonbar
        def ubar(g=gamma): return lu1*u1(g) + lu2*u2(g)
        def vbar(g=gamma): return ubar(g) - qbar(g)
        def ebar(g=gamma): return eps*cbar(g)/lc

        # Define rotated vectors
        def u(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@ubar(g)
        def e(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@ebar(g)
   
        # Derive angles for computing candidate length lw1 and
        # defining the edges correctly in Miura-ori based model
        beta1 = np.arccos(np.dot(-cbar(gamma0), vbar(gamma0))/\
                          (lc*np.linalg.norm(vbar(gamma0))))
        beta2 = np.arccos(np.dot(-cbar(gamma0), d(gamma0) - ls*a(gamma0))/\
                          (lc*np.linalg.norm(d(gamma0) - ls*a(gamma0))))
        beta3 = np.arccos(np.dot(-cbar(gamma0), b(gamma0) - ls*a(gamma0))/\
                              (lc*np.linalg.norm(b(gamma0) - ls*a(gamma0))))

        # Solve for additional vector lengths
        if ssa == 1:
            
            # Get the length of w
            lw1 = lc + lv1/np.tan(beta1) - lv1/np.tan(beta2)
            lw2 = lc - delta/np.tan(alpha1) - delta/np.tan(alpha2 - alpha1)
            lw = min(lw1, lw2)

        elif ssa == 2:

            # Track whether or not it's a symmetric case
            symmetric = np.allclose(alpha1 + alpha2, np.pi)
            
            if symmetric:
                
                # Redefine angles simply
                beta2 = alpha1
                beta3 = beta2

                # Solve for lengths of h and w
                lh = (lc - eps)/2
                lw = (lc - eps)/2

                # Get optimal length of wstar
                lwstar = lc - eps - delta/np.tan(alpha2 - alpha1)

            else:

                # Solve for important lengths using intersection
                p1 = c(gamma0)
                p2 = abar(gamma0) + cbar(gamma0) + dbar(gamma0) + s(gamma0) + ebar(gamma0) - ubar(gamma0)
                v1 = c(gamma0)/lc
                v2 = -cbar(gamma0)/lc

                lh = intersection(p1, p2, v1, v2)
                lw = intersection(p2, p1, v2, v1)

            # Get lengths of r and t
            lr = lu1/np.sin(beta3)
            lt = lv1/np.sin(beta2)

            # Define rbar and tbar and then solve for hbar
            def rbar(g=gamma): return -(lu1/np.tan(beta3) - lu2)*cbar(g)/lc + ubar(g)
            def tbar(g=gamma): return -(lv1/np.tan(beta2) + lv2)*cbar(g)/lc + vbar(g)
            def hbar(g=gamma): return lh*miuraori(rbar(g), cbar(g), tbar(g))

            # Define rotated vectors
            def h(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@hbar(g)
            def r(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@rbar(g)
            def t(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@tbar(g)

        # Define wbar and its rotation
        def wbar(g=gamma): return -lw*cbar(g)/lc
        def w(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@wbar(g)

        if ssa == 2:
            
            # Define lengths for f and g
            if symmetric:
                lfh = (2*lc - lwstar) - (lw + lu1/np.tan(beta3))
                lfr = 0
                lgh = lfh
                lgt = lfr

            #Asymmetric case
            else:
                #O1p and O3_n1
                p1 = abar(gamma0) + cbar(gamma0) + dbar(gamma0) + s(gamma0) + ebar(gamma0) - ubar(gamma0) + wbar(gamma0) + rbar(gamma0)
                p2 = b(gamma0) - d(gamma0) - e(gamma0) + u(gamma0) - w(gamma0)
                v1 = hbar(gamma0)/lh
                v2 = -rbar(gamma0)/lr
                            
                lfh = intersection(p1, p2, v1, v2)
                lfr = intersection(p2, p1, v2, v1)

                #Q1p and O_3
                p1 = abar(gamma0) + cbar(gamma0) + dbar(gamma0) + s(gamma0) + ebar(gamma0) - ubar(gamma0) + wbar(gamma0) + tbar(gamma0)
                p2 = -e(gamma0) + u(gamma0) - w(gamma0)
                v1 = hbar(gamma0)/lh
                v2 = -tbar(gamma0)/lt
                
                lgh = intersection(p1, p2, v1, v2)
                lgt = intersection(p2, p1, v2, v1)

            # Define fbar and gbar
            def fbar(g=gamma): return lfh*hbar(g)/lh + lfr*rbar(g)/lr
            def gbar(g=gamma): return lgh*hbar(g)/lh + lgt*tbar(g)/lt

            # Define rotated vectors
            def f(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@fbar(g)
            def g(g=gamma): return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])@gbar(g)

    """ DEFINE VERTICES AND EDGES """

    # Tube 1 Vertices
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

    # Tube 2 Vertices
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

    # Define all the edges in the zipper-coupled tubes
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

        # Points on adjacent tubes
        F0_1 = -b() + d()
        F3_1 = c() - b() + d()
        P0_1 = -e() - b() + d()
        
        F0p_n1 = s() - bbar() + dbar()
        F3p_n1 = s() + cbar() - bbar() + dbar()
        P0p_n1 = s() - ebar() - bbar() + dbar()

        if ssa == 1:
            
            # Smooth Sheet Attachment 1 Vertices
            H0 = u()
            H1 = a() + c() + d() - u() + w()
            H2 = a() + c() + d() - u()
            H3 = u() - w()

            # Smooth Sheet Attachment 2 Vertices
            H0p = s() + ubar()
            H1p = abar() + cbar() + dbar() + s() - ubar() + wbar()
            H2p = abar() + cbar() + dbar() + s() - ubar()
            H3p = s() + ubar() - wbar()

            # Define all the edges in the smooth sheet attachment
            sheet_edges = [(E1, H1), (E2, H2), (G1, H1), (G2, H2), # Sheet 1
                           (F0, H0), (F3, H3), (F0_1, H0), (F3_1, H3),
                           (H0, H3), (H1, H2),
                           
                           (E1p, H1p), (E2p, H2p), (G1p, H1p), (G2p, H2p), # Sheet 2
                           (F0p, H0p), (F3p, H3p), (F0p_n1, H0p), (F3p_n1, H3p),
                           (H0p, H3p), (H1p, H2p)]

        elif ssa == 2:
            
            # Smooth Sheet Attachment 1 Vertices
            O0 = -e() + u()
            O3 = -e() + u() - w()
            P0 = -e()

            O1 = a() + c() + d() + e() - u() + w() + r()
            O2 = a() + c() + d() + e()
            O4 = a() + c() + d() + e() - u() + w() + r() + f()

            P1 = a() + c() + d() + e() - u() + w()
            P2 = a() + c() + d() + e() - u()
            P4 = a() + c() + d() + e() - u() + w() + h()

            Q1 = a() + c() + d() + e() - u() + w() + t()
            Q2 = a() + b() + c() + e()
            Q4 = a() + c() + d() + e() - u() + w() + t() + g()

            # Smooth Sheet Attachment 2 Vertices
            O0p = s() - ebar() + ubar()
            O3p = s() - ebar() + ubar() - wbar()
            P0p = s() - ebar()

            O1p = abar() + cbar() + dbar() + s() + ebar() - ubar() + wbar() + rbar()
            O2p = abar() + cbar() + dbar() + s() + ebar()
            O4p = abar() + cbar() + dbar() + s() + ebar() - ubar() + wbar() + rbar() + fbar()

            P1p = abar() + cbar() + dbar() + s() + ebar() - ubar() + wbar()
            P2p = abar() + cbar() + dbar() + s() + ebar() - ubar()
            P4p = abar() + cbar() + dbar() + s() + ebar() - ubar() + wbar() + hbar()

            Q1p = abar() + cbar() + dbar() + s() + ebar() - ubar() + wbar() + tbar()
            Q2p = abar() + bbar() + cbar() + s() + ebar()
            Q4p = abar() + cbar() + dbar() + s() + ebar() - ubar() + wbar() + tbar() + gbar()

            # Define S1 and S3 differently if in the symmetric case
            if symmetric:
                O3 = -e() + u() - lwstar*w()/lw
                O3p = s() - ebar() + ubar() - lwstar*wbar()/lw

            # Define all the edges in the smooth sheet attachment
            sheet_edges = [(O0, P0), (O0, P0_1), (O3, F3), (O3, F3_1), # Sheet 1
                           (O0, O3), (P0, F3), (P0_1, F3_1),
                           (O1, P1), (O2, P2), (O4, P4),
                           (P1, Q1), (P2, Q2), (P4, Q4),
                           (O1, O2), (O1, O4), (P1, P2),
                           (P1, P4), (Q1, Q2), (Q1, Q4),
                           
                           (O0p, P0p), (O0p, P0p_n1), (O3p, F3p), (O3p, F3p_n1), # Sheet 2
                           (O0p, O3p), (P0p, F3p), (P0p_n1, F3p_n1),
                           (O1p, P1p), (O2p, P2p), (O4p, P4p),
                           (P1p, Q1p), (P2p, Q2p), (P4p, Q4p),
                           (O1p, O2p), (O1p, O4p), (P1p, P2p),
                           (P1p, P4p), (Q1p, Q2p), (Q1p, Q4p)]

        # Combine zipper-coupled tubes and smooth sheet attachment
        edges = np.concatenate((edges, sheet_edges))

    # Compose all the edges into a single list for plotting
    # Apply the transformation Phi
    component = np.empty((0, 6))
    for edge in edges:
        component = np.vstack((component, np.concatenate((Phi()@edge[0] + Phi()@origin,
                                                          Phi()@edge[1] - Phi()@edge[0]))))

    return component
