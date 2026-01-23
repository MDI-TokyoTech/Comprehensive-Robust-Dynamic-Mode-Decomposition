# Comprehensive Robust Dynamic Mode Decomposition from Mode Extraction to Dimnensional Reduction

Yuki Nakamura, Shingo Takemoto, and Shunsuke Ono

MDI Lab, Institute of Science, Tokyo, Japan


### Links
- [Preprint (arXiv)](https://arxiv.org/abs/2601.11116)

## Abstract
We propose Comprehensive Robust Dynamic Mode Decomposition (CR-DMD), a novel framework that robustifies the entire DMD process - from mode extraction to dimensional reduction - against mixed noise. Although standard DMD widely used for uncovering spatio-temporal patterns and constructing low-dimensional models of dynamical systems, it suffers from significant performance degradation under noise due to its reliance on least-squares estimation for computing the linear time evolution operator. Existing robust variants typically modify the least-squares formulation, but they remain unstable and fail to ensure faithful low-dimensional representations. First, we introduce a convex optimization-based preprocessing method designed to effectively remove mixed noise, achieving accurate and stable mode extraction. Second, we propose a new convex formulation for dimensional reduction that explicitly links the robustly extracted modes to the original noisy observations, constructing a faithful representation of the original data via a sparse weighted sum of the modes. Both stages are efficiently solved by a preconditioned primal-dual splitting method. Experiments on fluid dynamics datasets demonstrate that CR-DMD consistently outperforms state-of-the-art robust DMD methods in terms of mode accuracy and fidelity of low-dimensional representations under noisy conditions.

