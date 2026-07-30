# RxDAWN: Regularized xDAWN for ERP-Based Brain–Computer Interfaces

MATLAB implementation of the **Regularized xDAWN (RxDAWN)** algorithm proposed in:

> Mobaien A., Boostani R., Sanei S.  
> **Improving the performance of P300-based BCIs by mitigating the effects of stimuli-related evoked potentials through regularized spatial filtering.**  
> *Journal of Neural Engineering*, 21(1), 016023, 2024.  
> DOI: https://doi.org/10.1088/1741-2552/ad2495

---

## Overview

RxDAWN is a regularized extension of the classical xDAWN spatial filtering algorithm for enhancing P300 event-related potentials (ERPs) while simultaneously suppressing stimulus-related visual evoked potentials (VEPs). The method improves signal quality for ERP-based Brain–Computer Interface (BCI) applications by incorporating VEP suppression into the spatial filter optimization framework.

---

## Features

- MATLAB implementation of the RxDAWN algorithm
- Enhancement of P300 responses
- Suppression of stimulus-related VEPs
- Generalized eigenvalue formulation
- Compatible with ERP-based BCI research

---

## Requirements

- MATLAB R2018a or later (earlier versions may also work)
- Signal Processing Toolbox (if required by your workflow)

---

## Usage

```matlab
alpha = 0.5;
P300Length = 0.8;

[U, A_hat, eigenvalues] = RxDAWN( ...
    X, targetStimulusOnsets, Fs, P300Length, Cv, alpha);
```

---

## Inputs

| Parameter | Description |
|-----------|-------------|
| X | EEG data matrix (Nt × Ns) |
| targetStimulusOnsets | Sample indices of target stimuli |
| Fs | Sampling frequency (Hz) |
| P300Length | ERP duration (seconds) |
| Cv | VEP covariance matrix |
| alpha | Regularization parameter (0–1) |

---

## Outputs

| Output | Description |
|--------|-------------|
| U | Spatial filters |
| A_hat | Estimated ERP template |
| eigenvalues | Generalized eigenvalues |

---

## Citation

If you use this implementation, please cite:

```bibtex
@article{Mobaien2024,
  author = {Ali Mobaien and Reza Boostani and Saeid Sanei},
  title = {Improving the performance of P300-based BCIs by mitigating the effects of stimuli-related evoked potentials through regularized spatial filtering},
  journal = {Journal of Neural Engineering},
  volume = {21},
  number = {1},
  pages = {016023},
  year = {2024},
  doi = {10.1088/1741-2552/ad2495}
}
```

---

## License

Released under the MIT License.

---

## Author

Ali Mobaien  
