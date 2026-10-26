# The hybrid of configuration sets

The Discussion of the paper evaluates a hybrid of SWRT that applies one of three configuration sets to each spectrum, chosen by the class of the spectrum (the paper's Table 4). The hybrid was evaluated offline. NL6 reads one operator per run, and this version of the model does not change the operator during a run.

| Class | Rule | File | K | Offsets | S_nl relative to DIA |
|---|---|---|---|---|---|
| Opposing seas | r_w ≥ 5 % | `set_X.ww3op` | 33 | 19 | 3.7 |
| Narrow swell | otherwise σ_w < 20° | `set_S.ww3op` | 20 | 13 | 2.4 |
| Broad seas | all other spectra | `set_W.ww3op` | 33 | 18 | 3.6 |

For comparison, SWRT (`../swrt_K36.ww3op`, K = 36 on 18 offsets) takes 3.9 times as long as DIA.

Both quantities come from the directional distribution weighted by f⁴, D(θ) = Σ_f E(f, θ) f⁴ Δf, which stresses the high-frequency waves. σ_w is the circular directional spread of D, and r_w is the fraction of D more than 120° from the mean direction of D. `classify.py` gives the class and the file of each spectrum of a file. The limit of 20° is at the minimum of the distribution of σ_w over the synthetic spectra that are not opposing seas.

For each class, sets were selected by the search of SWRT (that of `tools/build_swrt.py` under an offset limit) on the synthetic training spectra of the class, along lines of equal modeled cost of S_nl. The set retained for a class is the cheapest with both a mean 2-D error ē2 and a median ε no larger than those of SWRT on the synthetic training spectra of the class.

The times of S_nl are those of Table 4 of the paper: W3SNL6 with each file against W3SNL1 (DIA), outside the model, on the 15,000 synthetic training spectra. The files have the format of `inputs/swrt_K36.ww3op` and belong to the same spectral grid, and each can be given to NL6 as its operator. With `tools/check_nl6.py`, W3SNL6 with each file reproduces the Python operator of its set to 8e-7 in S_nl and 4e-7 in the diagonal term.
