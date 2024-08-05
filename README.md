# TOV_Basics

This project provides a solution for understanding and solving the Tolman–Oppenheimer–Volkoff (TOV) equation, which describes the structure of a static, spherically symmetric body of isotropic material in general relativity. It includes Python scripts for handling different equations of state (EOS) and generating mass-radius relationships.

## Components

### Python Scripts:

- **PiecewiseEOS.py**:
  - Handles piecewise equations of state (EOS) for TOV calculations.
  - Uses interpolation methods (`interp1d`, `RegularGridInterpolator`, `LinearNDInterpolator`) for EOS data.
  - Provides a flexible way to define EOS with different segments.

- **PolytropicEOS.py**:
  - Handles polytropic equations of state for TOV calculations.
  - Suitable for modeling neutron stars and other astrophysical objects with a polytropic relation between pressure and density.

- **TabularEOS.py**:
  - Handles tabular equations of state for TOV calculations.
  - Reads EOS data from tables and uses interpolation to obtain values between table entries.

- **MR_Generator.py**:
  - Generates mass-radius relationships using different EOS models.
  - Utilizes the defined EOS scripts to solve the TOV equation and plot the mass-radius curve.

## Requirements

To run these scripts, you need to have the following Python packages installed:

- `numpy`
- `scipy`
- `matplotlib`
- `IPython`

You can install these packages using `pip`:

```bash
pip install numpy scipy matplotlib ipython
```
## Acknowledgments
Thanks to the developers of numpy, scipy, matplotlib, and IPython for their valuable libraries, which made this project possible.
