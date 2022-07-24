import numpy as np
import matplotlib.pyplot as plt
from AZCT import azct

def __getBounds(edges):
    """ Finds the x and y bounds for a given set of 2D edges

    Parameters:
        edges (numpy array): All the edges in a pattern--only four columns

    Returns:
        Tuple of bounds for plotting the structure
    """

    # Initialize empty sets for storing x and y coordinates
    Xs = set()
    Ys = set()

    # Loop through edges and get coordinates of both vertices in each edge
    for edge in edges:
        Xs.add(edge[0])
        Xs.add(edge[0] + edge[2])
        Ys.add(edge[1])
        Ys.add(edge[1] + edge[3])
        
    return min(Xs), max(Xs), min(Ys), max(Ys)


def printTube(edges, ssa, n, bounds):
    """ Prints origami folding pattern as a jpeg image
        for the asymmetric zipper-coupled tube given a set of edges

    Parameters:
        edges (numpy array): The edges in the origami structure
                             outputted by AZCT.azct
        n (int): The number of components in the structure
    """

    # Get the bounds for plotting
    xmin, xmax, ymin, ymax = bounds
    
    # Loop n times and derive
    mountain = np.empty((0,4))
    valley = np.empty((0,4))
    for i in range(n):

        # Get the basic cell and its complement from the vertices
        basic = edges[[40*i + x for x in [0, 3, 4, 7, 8, 11, 12,
                                          13, 15, 17, 16, 19]], :][:, [0, 2, 3, 5]]
        
        complement = edges[[40*i + x for x in [1, 2, 5, 6, 9, 10,
                                               13, 14, 15, 17, 18, 19]], :][:, [0, 2, 3, 5]]

        # Shift basic cells up and over and reflect
        for edge in basic:
            edge[1] += (.2 - ymin)
            edge[0] *= -1
            edge[2] *= -1
            edge[0] += (xmax + xmin)
            
        # Shift complementary cells down
        for edge in complement:
            edge[1] -= (.2 + ymax)

        # Define which vectors to keep
        if n == 1:
            keepM = [0, 1, 3, 4, 5, 6, 7, 8, 10, 11]
            keepV = [2]
            bold = np.vstack((basic[9], complement[9]))
        elif i == 0:
            keepM = [0, 3, 4, 5, 6, 7, 8, 10, 11]
            keepV = [1, 2]
            bold = np.vstack((basic[9], complement[9]))
        elif i < n - 1:
            keepM = [0, 3, 6, 7, 8, 9, 10, 11]
            keepV = [1, 2]
        else:
            keepM = [0, 1, 3, 6, 7, 8, 9, 10, 11]
            keepV = [2]

        # Add mountain and valley folds
        mountain = np.append(mountain,
                             np.concatenate((basic[keepM],
                                             complement[keepM])), axis=0)
        valley = np.append(valley,
                           np.concatenate((basic[keepV],
                                           complement[keepV])), axis=0)
    
    # Plot bold creases
    fig, ax = plt.subplots()
    X, Y, U, V = zip(*bold)
    ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1,
              linewidth=2, width=0.0001, headwidth=0, headlength=0,
              facecolor="none")

    # Plot mountain folds
    X, Y, U, V = zip(*mountain)
    ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1,
              linewidth=1, width=0.0001, headwidth=0, headlength=0,
              facecolor="none")

    # Plot valley folds
    X, Y, U, V = zip(*valley)
    ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1,
              linewidth=1, width=0.0001, headwidth=0, headlength=0,
              facecolor="none", linestyle=(0, (1, 3)))

    # Plot and save to file
    ax.set_xlim([xmin - .3, xmax + .3])
    ax.set_ylim([-(ymax - ymin + .5), ymax - ymin + .5])
    ax.axis('off')
    ax.set_aspect('equal')
    plt.savefig("AZCT Pattern.jpg", dpi=600)


def printSheet(edges, ssa, n, bounds):
    """ Prints origami folding pattern as a jpeg image
        for the smooth sheet attachment given a set of edges

    Parameters:
        edges (numpy array): The edges in the origami structure
                             outputted by AZCT.azct
        n (int): The number of components in the structure
    """

    # Get the bounds for plotting
    xmin, xmax, ymin, ymax = bounds

    # Specify the indices of the edges to be plotted in cell1 and cell2
    # Also specify which edges are valley folds and which are mountain
    if ssa == 1:
        num = 60
        
        indices1 = [21, 29, 50, 51, 52, 53, 59]
        m1 = [0, 1, 2, 3, 4, 5]
        v1 = [6]
        
        indices2 = indices1
        m2 = m1
        v2 = v1
        
    elif ssa == 2:
        num = 78
        shift = edges[4][3]
        
        indices1 = [66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77]
        m1 = [1, 2, 4, 5, 6, 7, 9, 10, 11]
        v1 = [0, 3, 8]

        indices2 = [40, 41, 42, 43, 44, 45, 46]
        m2 = [0, 1, 2, 3, 5, 6]
        v2 = [4]

    # Loop n times and derive 
    mountain = np.empty((0,4))
    valley = np.empty((0,4))
    for i in range(n):

        # Get the cell edges from the structure
        cell1 = edges[[num*i + x for x in indices1]][:, [0, 1, 3, 4]]

        # Define second cell which is a reflection of the first
        cell2 = edges[[num*i + x for x in indices2]][:, [0, 1, 3, 4]]

        # Flip cell2 if regular and then shift all cells
        if ssa == 1:
            for edge in cell2:
                edge[1] *= -1
                edge[3] *= -1
        
            for edge in cell1:
                edge[1] -= (ymax + .2)
            for edge in cell2:
                edge[1] += (ymax + .2)

        elif ssa == 2:
            for edge in cell1:
                edge[1] -= (ymax)
            for edge in cell2:
                edge[1] -= (ymax - .4)
            
        # Add mountain and valley folds
        mountain = np.append(mountain, np.concatenate((cell1[m1],
                                             cell2[m2])), axis=0)
        valley = np.append(valley, np.vstack((cell1[v1],
                                             cell2[v2])), axis=0)
    
    # Plot mountain folds
    fig, ax = plt.subplots()
    X, Y, U, V = zip(*mountain)
    ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1,
              linewidth=1, width=0.0001, headwidth=0, headlength=0,
              facecolor="none")

    # Plot valley folds
    X, Y, U, V = zip(*valley)
    ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1,
              linewidth=1, width=0.0001, headwidth=0, headlength=0,
              facecolor="none", linestyle=(0, (1, 3)))

    # Plot and save to file
    ax.set_xlim([xmin - .3, xmax + .3])
    ax.set_ylim([-(ymax - ymin + .5), (ymax - ymin + .5)])
    ax.axis('off')
    ax.set_aspect('equal')
    plt.savefig("SSA Pattern.jpg", dpi=600)
    

# Main loop for user interface
if __name__ == '__main__':

    # Number of components
    n = int(input("Enter the number of components in Z,\n" +
                  "your asymmetric zipper-coupled tubes structure: "))
    if n < 1:
        raise ValueError("Number must be a positive integer")

    # Whether or not to have non-uniform extensions
    nue = 'n'
    if n > 1:
        nue = input("Does your structure contain non-uniform extensions? (y/n): ")
        if nue not in ['y', 'n']:
            raise ValueError("Please respond either y or n")
            
    # Whether or not to have smooth sheet attachments and which type
    ssa = 0
    if nue != 'y':
        if input("Does your structure have smooth sheet attachments? (y/n): ") == 'y':
            ssa = int(input("Are the attachments regular (1) or Miura-ori inspired (2)?: "))
            if ssa not in [1, 2]:
                raise ValueError("Attachment type must be either 1 or 2")
            
    origin_flat = np.array([0, 0, 0])
    origin_ideal = np.array([0, 0, 0])
    lb0 = None

    # Loop through each component
    for i in range(n):

        # Get length of b if non-uniform extension
        if nue == 'y' and i > 0:
            if i == 1:
                lb0 = lb
            lb = float(input("  Length of crease b in component "
                        + str(i + 1) + ": "))
            if lb <= 0:
                raise ValueError("Lengths must be positive")
            
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
            
            # Set gamma and find gamma0
            gamma = 1e-5
            if ssa:
                ax = (np.sin(alpha1)*np.cos(alpha4) - np.sin(alpha4)*np.cos(alpha1))
                bx = (np.sin(alpha1)*np.cos(alpha3) - np.sin(alpha4)*np.cos(alpha2))
                ay = (np.sin(alpha1)*np.cos(alpha4) + np.sin(alpha4)*np.cos(alpha1))
                by = (np.sin(alpha1)*np.cos(alpha3) + np.sin(alpha4)*np.cos(alpha2))
                k = 2*np.sin(alpha1)*np.sin(alpha4)
                gamma0 = np.arccos(np.sqrt((ay*by - ax**2 + k**2 + \
                                            np.sqrt((ay*by + ax**2 + k**2)**2 - \
                                                    4*(ay**2)*ax*bx))/(2*(k**2)))) + 1e-5

                print("\nNote: if a runtime warning was thrown, invalid values" +
                      "\nfor alpha1 and alpha2 were used--refer to Figure 11\n")

                # Set bhat and dhat for ideal state amd get the edges
                bhat_ideal = np.array([-np.sin(alpha1)*np.cos(gamma0),
                                      np.sin(alpha1)*np.sin(gamma0),
                                      -np.cos(alpha1)])
                dhat_ideal = np.array([np.sin(alpha4)*np.cos(gamma0),
                                      np.sin(alpha4)*np.sin(gamma0),
                                      -np.cos(alpha4)])
                edges_ideal = azct(origin_ideal, alpha1, alpha2, la, lb, lc, gamma0, ssa=ssa)

            
            # Set bhat and dhat for unfolded state
            bhat_flat = np.array([-np.sin(alpha1)*np.cos(gamma),
                                  np.sin(alpha1)*np.sin(gamma),
                                  -np.cos(alpha1)])
            dhat_flat = np.array([np.sin(alpha4)*np.cos(gamma),
                                  np.sin(alpha4)*np.sin(gamma),
                                  -np.cos(alpha4)])

            # Define the first component
            edges_flat = azct(origin_flat, alpha1, alpha2, la, lb, lc, gamma, ssa=0, lb0=lb0)
            
        else:
            # Update origin and attach new component in unfolded state
            origin_flat = origin_flat + ld*dhat_flat - lb*bhat_flat
            component = azct(origin_flat, alpha1, alpha2, la, lb, lc, gamma, ssa=0, lb0=lb0)
            edges_flat = np.concatenate((edges_flat, component))

            if ssa:
                # Update origin and attach new component in ideal state
                origin_ideal = origin_ideal + ld*dhat_ideal - lb*bhat_ideal
                component = azct(origin_ideal, alpha1, alpha2, la, lb, lc, gamma0, ssa=ssa)
                edges_ideal = np.concatenate((edges_ideal, component))

        # Update the length of d in the most recent component
        ld = lb*(np.sin(alpha1)/np.sin(alpha4))

    # Get the bounds for plotting
    xmin1, xmax1, ymin1, ymax1 = __getBounds(edges_flat[:, [0, 2, 3, 5]])
    xmin2, xmax2, ymin2, ymax2 = __getBounds(edges_ideal[:, [0, 1, 3, 4]])
    bounds = (min(xmin1, xmin2), max(xmax1, xmax2), min(ymin1, ymin2), max(ymax1, ymax2))
    
    # Plot the resulting structure
    printTube(edges_flat, ssa, n, bounds)
    if ssa:
        printSheet(edges_ideal, ssa, n, bounds)
