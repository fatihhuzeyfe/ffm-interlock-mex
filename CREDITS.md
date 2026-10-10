# Credits

## Material data

* Elastic constants and strengths of PLA, PETG and ABS: UltiMaker technical data sheets
  (v5.00, test protocol ASTM D3039 in XY, YZ and Z orientation).
* Orientation-dependent toughness of printed PLA: Kiendl and Gao, Compos. Part B 180 (2020) 107562.
* Fracture toughness of printed ABS: Hart and Wetzel, Eng. Fract. Mech. 177 (2017) 1-13.
* Details and assumptions: [docs/MATERIALS.md](docs/MATERIALS.md).

## Methods

The coupled criterion follows Leguillon (2002) and Cornetti et al. (2006); the mixed-mode
partition follows the crack-closure approach of Rybicki and Kanninen (1977); mesh convergence is
reported with the grid convergence index of Roache (1994) and Celik et al. (2008).

## Use of AI

The code, the figures and the documentation of this repository were developed with the
assistance of Claude (Anthropic). All finite element analyses were run and all results were
verified by the author.
