# ARBIFOIL
Airfoil analysis library

_Theodorsen's potential theory of arbitrary wing sections implemented in Python 3.7._

## Using Arbifoil
Generate or otherwise obtain an airfoil .dat file in **Selig format***.

Import arbifoil.py 
```
from arbifoil import *
```

Instantiate a foil object by passing it a .dat airfoil coordinates file. 
```
clarky = foil('airfoils/clarky.txt')
```

After the mapping is determined, you have at your disposal the following methods for analysis:
- **C_L(aoa)** - returns the coefficient of lift at specified angle of attack in degrees,
- **C_M_AC()** - returns the moment coefficient about aerodynamic center,
- **C_M_LE(aoa)** - returns the moment coefficient about leading edge at specified angle of attack in degrees,
- **C_M_1Q(aoa)** - returns the moment coefficient about 1 quarter chord point at specified angle of attack in degrees,
- **CoP(aoa)** - returns nondimensional center of pressure measured from the leading edge,
- **AC()** - returns the nondimensional coordinates of the aerodynamic center.

## DAT File Rules

**Selig format** 

The points run from the trailing edge across the upper surface of the airfoil
to the leading edge and back across the lower surface. 
For examples, check the arbifoil/airfoils/ folders.

**Coordinate system positioning**

The origin of the _right handed_ system is at the leading edge. 

The x axis points towards the trailing edge.

The upper surface points have positive y values, the lower surface points have negative y values.

## References
[1] Theodorsen, T.: Theory of Wing Sections of Arbitrary Shape, NACA, 1931.

[2] Theodorsen, T., Garrick I.E.: General Potential Theory of Arbitrary Wing Sections, NACA TR-452, 1934.

[3] Karamcheti, K.: Principles of ideal-fluid aerodynamics, R. E. Krieger Publishing Company, 1980.

