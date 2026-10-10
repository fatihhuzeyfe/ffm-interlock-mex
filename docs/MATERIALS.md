# Material data

Axes: 1 = road direction, 2 = in-plane transverse (road-to-road bond), 3 = build direction.
The 2D model lies in the layer plane (1-2), plane strain in direction 3. Units: MPa, N/mm (= kJ/m^2).

## Bulk materials of the study

| Material | E1 | E2 = E3 | nu12 | nu23 | G12 | G23 | X | Y | S | G_L | G_T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PLA | 3292 | 3071 | 0.35 | 0.35 | 1188.1 | 1137.4 | 56.0 | 33.1 | 27.8 | 1.5 | 4.5 |
| PETG | 1874 | 1711 | 0.38 | 0.38 | 656.8 | 619.9 | 50.3 | 19.0 | 16.0 | 1.5 | 8.0 |
| ABS | 1931 | 1699 | 0.35 | 0.35 | 681.8 | 629.3 | 38.2 | 19.0 | 16.0 | 0.5 | 2.26 |

X, Y: tensile strength along / across the roads; S: in-plane shear strength; G_L, G_T: toughness
for cracks across / along the roads. Sources: UltiMaker technical data sheets (elastic constants,
strengths); Kiendl and Gao (2020) for PLA and Hart and Wetzel (2017) for ABS (toughness); PETG
toughness estimated (least certain input).

## Reference interfaces (k_sigma = k_G = 1)

Rule `weaker_intralayer`: strength and toughness of the road-to-road bond of the weaker material;
sigma_c = 19 MPa, tau_c = 16 MPa; G_Ic = 1.5 kJ/m^2 (PLA/PETG), 0.5 kJ/m^2 (pairs with ABS);
G_IIc / G_Ic in {1, 2, 4}. In the study, strength and toughness are scaled by k_sigma and k_G.
