import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from AZCT import azct

def plot(edges, ssa):
    """ Plots edges using matplotlib in 3D

    Parameters:
        edges (2D numpy array): Contains the edges to plot (rows),
                formatted as vectors of length six: the first three
                components are tail position and the next three direction
        ssa (string): String signifying whether or not there are edges
                in a smooth sheet attachment to color
    """
    
    # Define colors
    cols = np.tile(np.array([0, 0, 0]), (len(edges), 1))
    if ssa:
        for i in range(len(edges)//60):
            cols[(60*i + 40):(60*i + 60)] = np.array([1, 0, 0])

    # Plot vectors
    X, Y, Z, U, V, W = zip(*edges)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(X, Y, Z, U, V, W, colors=cols, arrow_length_ratio=0)

    # Set axes limits
    ax.set_xlim([np.amin(edges[:,0]) - 1, np.amax(edges[:,0]) + 1])
    ax.set_ylim([np.amin(edges[:,1]) - 1, np.amax(edges[:,1]) + 1])
    ax.set_zlim([np.amin(edges[:,2]) - 2, np.amax(edges[:,2]) + 2])
    #ax.axis("off")

    plt.show()

# Main loop for user interface
if __name__ == '__main__':

    n = int(input("Enter the number of components in Z,\n" +
                  "your asymmetric zipper-coupled tubes structure: "))

    # Whether or not to have non-uniform extensions
    nue = 'n'
    if n > 1:
        nue = input("Does your structure contain non-uniform extensions? (y/n): ")

    # Whether or not to have smooth sheet attachments
    ssa = False
    if nue != 'y':
        ssa = input("Does you structure have smooth sheet attachments? (y/n): ") == 'y'

    origin = np.array([0, 0, 0])
    lb0 = None

    # Loop through each component
    for i in range(n):

        # Get length of b if non-uniform extension
        if nue == 'y' and i > 0:
            if i == 1:
                lb0 = lb
            lb = float(input("  Length of crease b in component "
                        + str(i + 1) + ": "))

        # Get design angles and lengths on first loop
        if i == 0:

            # Get design angles and lengths as input
            print("\nPlease enter the following values (angles in degrees)")
            alpha1 = float(input("  Design angle alpha1: "))*np.pi/180
            alpha2 = float(input("  Design angle alpha2: "))*np.pi/180
            la = float(input("  Length of crease a: "))
            lb = float(input("  Length of crease b: "))
            lc = float(input("  Length of crease c: "))

            # Check for invalid entries
            if alpha1 > alpha2:
                raise ValueError("alpha1 must be less than alpha2")
            if alpha1 + alpha2 > np.pi:
                raise ValueError("alpha1 + alpha2 must be less than 180")
            if alpha1 <= 0 or alpha2 <= 0 or la <= 0 or lb <= 0 or lc <= 0:
                raise ValueError("angles and lengths must be positive")

            # Set alpha3 and alpha4
            alpha3 = np.pi - alpha1
            alpha4 = np.pi - alpha2
            
            # Find gamma0
            gamma0 = 0
            if ssa:
                ax = (np.sin(alpha1)*np.cos(alpha4) - np.sin(alpha4)*np.cos(alpha1))
                bx = (np.sin(alpha1)*np.cos(alpha3) - np.sin(alpha4)*np.cos(alpha2))
                ay = (np.sin(alpha1)*np.cos(alpha4) + np.sin(alpha4)*np.cos(alpha1))
                by = (np.sin(alpha1)*np.cos(alpha3) + np.sin(alpha4)*np.cos(alpha2))
                k = 2*np.sin(alpha1)*np.sin(alpha4)
                gamma0 = np.arccos(np.sqrt((ay*by - ax**2 + k**2 + \
                                            np.sqrt((ay*by + ax**2 + k**2)**2 - \
                                                    4*(ay**2)*ax*bx))/(2*(k**2))))*180/np.pi

                print("\nNote: if a runtime warning was thrown, invalid values" +
                      "\nfor alpha1 and alpha2 were used--refer to Figure 11\n")

            # Get gamma as input
            gamma = float(input("  Driving angle gamma (between "
                                + str(round(gamma0 + 1e-4, 4)) +
                                " and 90: "))*np.pi/180 - 1e-6

            # Check if gamma is valid
            if gamma*180/np.pi < gamma0 or gamma*180/np.pi > 90:
                raise ValueError("Incorrect value entered for gamma")
            
            # Set bhat and dhat
            bhat = np.array([-np.sin(alpha1)*np.cos(gamma),
                             np.sin(alpha1)*np.sin(gamma),
                             -np.cos(alpha1)])
            dhat = np.array([np.sin(alpha4)*np.cos(gamma),
                             np.sin(alpha4)*np.sin(gamma),
                             -np.cos(alpha4)])

            # Define the first component
            edges = azct(origin, alpha1, alpha2, la, lb, lc, gamma, ssa, lb0)

        else:
            # Update origin and attach new component
            origin = origin + ld*dhat - lb*bhat
            component = azct(origin, alpha1, alpha2, la, lb, lc, gamma, ssa, lb0)
            edges = np.concatenate((edges, component))

        # Update the length of d in the most recent component
        ld = lb*(np.sin(alpha1)/np.sin(alpha4))

    # Plot the resulting structure
    plot(edges, ssa)
