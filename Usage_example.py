"""This script demonstrates some very basic usage of the
Inertial Lift Force Helper Classes associated with the
Curved Trapezoidal Duct data (or ILFHC-CTD for short).
The acronyms VST (vertically symmetric trapezoid) and
FBT (flat-bottomed trapezoid) are used to distinguish
the two sets of classes/data.
The class methods are themselves documented via
```help(InertialLiftForceHelperVST)``` and
```help(InertialLiftForceHelperFBT)```.
Some effort has been made to ensure the interfaces are more
or less identical, while the inner workings may differ.

Note: use of these helper classes requires the following files
to be located relative to the working directory:

 - 'VST_data/strap_2x1_lift_data.npz'
 - 'VST_data/strap_4x1_lift_data.npz'
 - 'FBT_data/utrap_4x1_lift_data.npz'

Observe that data for two aspect ratios W/H=2,4 is
available for VSTs, but only W/H=4 is available for FBTs.

Please cite our research if you use this code/data.
The model is developed in a JFM paper
(doi.org/10.1017/jfm.2019.323).
An in-depth study of the dynamics within trapezoidal
cross-sections is explored in (TODO: add details).

This code is provided under an MIT license
(see https://opensource.org/licenses/MIT).
However, I ask that the associated research data not be
directly distributed to third parties, rather I would
appreciate if it was cloned directly from my own respository.
Please don't hesitate to contact me if you have any queries.

Brendan Harding, 2023."""

# Import the helper class for vertically symmetric trapezoid cross-sections
from ILFHC_VST import InertialLiftForceHelperVST

# Initialise an instance for a neutrally buoyant particle with radius a=0.10
# within a curved duct having bend radius R=80.0 and a vertically symmetric
# trapezoid cross-section with aspect ratio 2 and shape parameter D=0.2.
VSTH = InertialLiftForceHelperVST(0.10,80.0,2,0.2)

# Print an estimate of the migration force at the
# centre of the cross-section (includes both the
# inertial lift and secondary flow contributions)
print(VSTH.migration_force(0.0,0.0))

# Generate a plot of the migration force
VSTH.plot_migration_force()

# Import the helper class for flat-bottomed trapezoid cross-sections
from ILFHC_FBT import InertialLiftForceHelperFBT

# Initialise an instance for a neutrally buoyant particle with radius a=0.10
# within a curved duct having bend radius R=160.0 and a flat-bottomed
# trapezoid cross-section with aspect ratio 4 and shape parameter D=0.2.
FBTH = InertialLiftForceHelperFBT(0.10,160.0,4,0.2)

# Print an estimate of the migration force at the
# centre of the cross-section (includes both the
# inertial lift and secondary flow contributions)
print(FBTH.migration_force(0.0,0.0))

# Generate a plot of the migration force
FBTH.plot_migration_force()

# Print the documentation for the class
# (including a description of the methods/functions available)
# (similar applies to the VST class)
help(InertialLiftForceHelperFBT) # or help(FBTH)
