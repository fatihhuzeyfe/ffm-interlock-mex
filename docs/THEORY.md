# Theory in brief

## Coupled stress-energy criterion

A crack of length Δa forms along a path Γ when two conditions hold simultaneously:

* stress condition, at every point of the future crack (interface: quadratic criterion)

  (⟨σ_n⟩/σ_c)² + (τ/τ_c)² ≥ 1   for all x in Γ(Δa);

  for bulk paths the normal stress on the crack plane is compared with the
  orientation-dependent strength of the printed material;

* energy condition, with the incremental energy release rate Ḡ = −ΔΠ/Δa and a linear
  mixed-mode law

  Ḡ_I/G_Ic + Ḡ_II/G_IIc ≥ 1.

Because the problem is linear, one FE solution per crack length suffices: stresses scale with
the load factor λ and energies with λ². With λ_s(Δa) and λ_G(Δa) the load factors that satisfy the
two conditions, the onset load of a path is

  λ_c = min over Δa of max(λ_s(Δa), λ_G(Δa)),

and the predicted onset is the path with the smallest λ_c. The implementation is in
`ffm/kernel.py` (`ffm_path`); Ḡ comes from strain-energy differences, its mode partition from
crack closure with the nodal forces of the uncracked model (`ffm/postprocess.py`).

## Interface length scale

With all lengths scaled by w, the stresses are self-similar and Ḡ = σ² w ĝ(Δa/w) / E*. The onset
stress relative to the interface strength then depends on the materials and the size only through

  l_ch / w,   l_ch = E* G_Ic / σ_c²,

besides the geometry, the elastic mismatch, τ_c/σ_c and G_IIc/G_Ic. Interface quality factors
enter as l_ch = (k_G / k_σ²) l_ch,ref.

* l_ch/w → 0: the stress condition is critical; at a re-entrant corner with singularity order
  λ − 1 the onset load decreases as (l_ch/w)^(1−λ), so sharp interlocks can be weaker than a
  flat interface.
* l_ch/w → ∞: the energy condition is critical and the crack length is limited by the interface
  segment between two features; the gain over the flat interface tends to a constant that depends
  only on the geometry (plateau of the master curve).

## Cohesive-zone counterpart

`scripts/make_czm_cases.py` builds the same unit cells with a cohesive interface (quadratic
traction initiation with σ_c, τ_c; linear softening with G_Ic, G_IIc; power law with exponent 1).
The CZM onset is the 5 % deviation from the initial linear response. When the cohesive length
approaches the size of the interface segment, the CZM predicts lower onset loads than the FFM;
the FFM gain is then an upper bound.

References: Leguillon (2002) Eur J Mech A/Solids 21:61; Cornetti et al. (2006) Eng Fract Mech
73:2021; Weißgraeber, Leguillon, Becker (2016) Arch Appl Mech 86:375.
