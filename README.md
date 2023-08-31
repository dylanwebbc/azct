# azct
Code for visualizing a 3D model of, and printing origami patterns for, asymmetric zipper-coupled tubes and their smooth sheet attachment. Below is a brief description of each python file. 

Read the academic journal articles here:
1. https://doi.org/10.1115/1.4056868
2. https://doi.org/10.3390/math10152643
3. https://doi.org/10.1115/DETC2022-90045 (conference publication of 1)

### AZCT
  This file contains the mathematics for modelling an asymmetric zipper-coupled tubes structure and its smooth sheet attachment. It does not have an interface, but may be of interest to the more serious user for its example of how to program the design process for the device.
  
### model
  This file allows the user to visualize a 3D model of the device---simply open the file and a UI will walk the user through the various necessary parameters.
  
### patterns
  This file allows the user to print out origami folding patterns for the device---simply open the file and a UI will walk the user through the various necessary parameters. There are two output files. AZCT Pattern.jpeg contains the patterns for a single asymmetric origami tube (the bolded line indicates where the patterns should be connected); print two copies of these patterns to have two tubes to couple together. SSA Pattern.jpeg contains the patters for the smooth sheet attachment for a single asymmetric origami tube; again, print two copies of these patterns to have enough smooth sheets for the structure.
