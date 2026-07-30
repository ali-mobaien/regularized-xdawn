# Regularized xDAWN

MATLAB implementation of the **Regularized xDAWN** spatial filtering algorithm proposed in:

> Mobaien A., Boostani R., Sanei S.
> *Improving the performance of P300-based BCIs by mitigating the effects of stimuli-related evoked potentials through regularized spatial filtering.*
> Journal of Neural Engineering, 21(1), 016023, 2024.
> https://doi.org/10.1088/1741-2552/ad2495

## Overview

Regularized xDAWN is a spatial filtering method designed to improve the detection of P300 event-related potentials (ERPs) in EEG recordings. The method extends the classical xDAWN algorithm by incorporating regularization to simultaneously enhance P300 responses while suppressing stimulus-related visual evoked potentials (VEPs), leading to improved performance in ERP-based brain–computer interfaces.

## Features

- MATLAB implementation
- Regularized spatial filtering
- P300 enhancement
- VEP suppression
- ERP-based BCI applications

## Citation

If you use this code in your research, please cite:

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

## License

MIT License.

## Author

Ali Mobaien
