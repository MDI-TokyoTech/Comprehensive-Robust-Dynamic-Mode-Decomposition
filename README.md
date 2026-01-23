# Comprehensive Robust Dynamic Mode Decomposition from Mode Extraction to Dimnensional Reduction

Yuki Nakamura, Shingo Takemoto, and Shunsuke Ono

MDI Lab, Institute of Science, Tokyo, Japan


### Links
- [Preprint (arXiv)](https://arxiv.org/abs/2601.11116)

## Abstract
We propose Comprehensive Robust Dynamic Mode Decomposition (CR-DMD), a novel framework that robustifies the entire DMD process - from mode extraction to dimensional reduction - against mixed noise. Although standard DMD widely used for uncovering spatio-temporal patterns and constructing low-dimensional models of dynamical systems, it suffers from significant performance degradation under noise due to its reliance on least-squares estimation for computing the linear time evolution operator. Existing robust variants typically modify the least-squares formulation, but they remain unstable and fail to ensure faithful low-dimensional representations. First, we introduce a convex optimization-based preprocessing method designed to effectively remove mixed noise, achieving accurate and stable mode extraction. Second, we propose a new convex formulation for dimensional reduction that explicitly links the robustly extracted modes to the original noisy observations, constructing a faithful representation of the original data via a sparse weighted sum of the modes. Both stages are efficiently solved by a preconditioned primal-dual splitting method. Experiments on fluid dynamics datasets demonstrate that CR-DMD consistently outperforms state-of-the-art robust DMD methods in terms of mode accuracy and fidelity of low-dimensional representations under noisy conditions.

## How to use

### Quick Start
The easiest way to use the code is to run `examples/main_cylinder.m`.
This script handles path setup, data loading, adding synthetic noise, running CR-DMD, and visualization all at once.

```matlab
% Run in MATLAB Command Window
>> run examples/main_cylinder.m
```

### Configuration (Parameters)
The `CR_DMD` class is configured via a `params` structure.
`params` consists mainly of three sections: (`pre`, `me`, `dr`).

```matlab
params.pre  % Parameters for Preprocessing
params.me   % Parameters for Mode Extraction
params.dr   % Parameters for Dimensional Reduction
crdmd = CR_DMD(params);
```

#### Preprocessing Parameters (`params.pre`)
| Parameter | Description | Default |
| :--- | :--- | :--- |
| `w` | Weight parameter for spatial differentiation | 0.5 |
| `epsilon` | Constraint parameter for noise level | 0.1 |
| `eta` | Constraint parameter for outlier sparsity (L1 norm) | 0.1 |
| `max_iter` | Maximum number of iterations | 5000 |
| `criteria` | Convergence criteria | 1e-4 |

#### Mode Extraction Parameters (`params.me`)
| Parameter | Description | Default |
| :--- | :--- | :--- |
| `r` | Target rank (number of modes to extract) | 21 |

#### Dimensional Reduction Parameters (`params.dr`)
| Parameter | Description | Default |
| :--- | :--- | :--- |
| `w` | Weight parameter for spatial differentiation | 0.5 |
| `epsilon` | Constraint parameter for noise level | 0.1 |
| `eta` | Constraint parameter for outlier sparsity | 0.1 |
| `mu` | Weight for mode selection sparsity | 1 |
| `max_iter` | Maximum number of iterations | 5000 |
| `criteria` | Convergence criteria | 1e-4 |

### Usage Details

#### 1. Initialization
```matlab
crdmd = CR_DMD(params);
```

#### 2. Preprocessing & Mode Extraction
Execute preprocessing to remove noise and extract DMD modes.
```matlab
[X_clean, S_outliers, Phi, Lambda, b] = crdmd.preprocess(X_noisy, h, v);
```
- **Inputs**:
    - `X_noisy`: Noisy data matrix ($N \times M$).
    - `h`, `v`: Spatial dimensions ($N = h \times v$).
- **Outputs**:
    - `X_clean`: Data after noise removal.
    - `S_outliers`: Detected sparse outliers.
    - `Phi`: DMD modes.
    - `Lambda`: Eigenvalues.
    - `b`: Mode amplitudes.

#### 3. Dimensional Reduction
Execute the dimensional reduction step to select the most significant modes.
```matlab
[xi, S_outliers_dr] = crdmd.dimensional_reduction(X_noisy, Phi, Lambda, b, h, v);
```
- **Inputs**:
    - `X_noisy`, `Phi`, `Lambda`, `b`: Results from the previous step.
    - `h`, `v`: Spatial dimensions.
- **Outputs**:
    - `xi`: Sparse coefficient vector selecting significant modes.
    - `S_outliers_dr`: Outliers detected in this stage.


