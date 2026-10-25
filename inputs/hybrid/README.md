# The hybrid of configuration sets

The Discussion of the paper evaluates a hybrid of SWRT that applies one of three configuration sets to each spectrum, chosen by the class of the spectrum (the paper's Table 4). The hybrid was evaluated offline. NL6 reads one operator per run, and this version of the model does not change the operator during a run.

| Class | Rule | File | K | Offsets |
|---|---|---|---|---|
| Opposing seas | r ≥ 5 % | `set_X.ww3op` | 28 | 16 |
| Narrow swell | otherwise σ_w < 20° | `set_S.ww3op` | 20 | 13 |
| Broad seas | all other spectra | `set_W.ww3op` | 33 | 20 |

σ_w is the circular directional spread of D(θ) = Σ_f E(f, θ) f⁴ Δf, and r is the fraction of the energy more than 120° from the mean direction of D (`classify.py`, which gives the class of each spectrum of a file). The f⁴ weight stresses the high-frequency waves. The limit of 20° is at the minimum of the distribution of σ_w over the synthetic spectra that are not opposing seas.

For each class, sets were selected by the search of SWRT (that of `tools/build_swrt.py` under an offset limit) on the synthetic spectra of the class, along lines of equal modeled cost of S_nl. The set retained for a class is the cheapest whose error on the spectra of the class in two new draws of synthetic spectra is no larger than that of SWRT. The files have the format of `inputs/swrt_K36.ww3op` and belong to the same spectral grid, and each can be given to NL6 as its operator.
