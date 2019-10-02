# arbifoil
Airfoil analysis library 
_Theodorsen's potential theory of arbitrary wing sections implemented in Python 3.7._

## Using Arbifoil
Generate or otherwise obtain a .dat point coordinates file in **Selig format***

Import arbifoil.py 
```
from arbifoil import *
```

Instantiate a foil object by passing it a .dat airfoil coordinates file. 
```
clarky = foil('airfoils/clarky.txt')
```

After the mapping is determined, you have at your disposal the following methods for analysis:
- C_L(aoa) - returns the coefficient of lift at specified angle of attack in degrees,
- C_M_AC() - returns the moment coefficient about aerodynamic center,
- C_M_LE(aoa) - returns the moment coefficient about leading edge at specified angle of attack in degrees,
- C_M_1Q(aoa) - returns the moment coefficient about 1 quarter chord point at specified angle of attack in degrees,
- CoP(aoa) - returns center of pressure nondimensionalized with the chord measured from the leading edge


***IMPORTANT***: Airfoil file has to be in **Selig format** with the points running from the trailing edge across
the upper surface of the airfoil to the leading edge and back across the lower surface.

