Linearized Supersonic Potential Flow Wing Solver

C# Cmdline utility that implements part of the methods of NASA TP 1718 to solve the flow over a supersonic wing, assuming linearized, irrotational flow.  
Intended to act as a proof of concept for implementing these methods into FAR at a later date.

Wing planforms can be defined using text files, using the format:

```
TRIANGLE
x1,y1,z1
x2,y2,z2
x3,y3,z3
```

Additional spaces or empty lines will cause errors, as will typos.  
The wing is assumed to be primarily in the Z plane; FUTURE IMPLEMENTATION: use Z values to account for wing camber and dihedral

Sideslip and AoA can be simulated as well

Simulation will dump a sim.csv file including the local Cp differences over the upper and lower surfaces that can be used to see the lift distribution  
The wing shape will be distorted as necessary to make the diagonals in the resulting spreadsheet mach lines