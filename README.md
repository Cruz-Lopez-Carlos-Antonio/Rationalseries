# Rational Series Resolutions Applied to Bateman Equations

## Overview of the Repository

The present repository contains the **Python 3** and **Wolfram Mathematica** codes associated with the evaluation of the analytical closed-relationship developed in our research for rational series of the exact form:

$$\sum_{k=0}^\infty \frac{P(k)\, z^k}{(k+a_1)^{m_1+1} (k+a_2)^{m_2+1} \cdots (k+a_n)^{m_n+1}}.$$

The framework extends beyond the explicit evaluation of these series involving rational terms; it constitutes a computable implementation of confluent divided differences (divided differences with repeated arguments). The set of codes accompany the manuscript *On the Generalized Summation of Series with Rational Coefficients*, recently submitted to the *Computer Physics Communications* journal. Unless otherwise noted, all scripts are released under the **MIT License**.

**Authors:**
* Carlos-Antonio Cruz-López
* Marc Jornet
* Gilberto Espinosa-Paredes
* Juan-Luis François

## Physical Application
A direct and strong physical application emerges by adapting this mathematical structure to solve the generalized Bateman equations. By mapping the algebraic-combinatorial structure of the rational series resolution to the general solution of the Bateman equations, the algorithm can directly and systematically compute complex decay chain and transmutation models, including cases with repeated decay constants.

## Theoretical Framework

The core mathematical foundation of this computational implementation is established by **Theorem 9** of our manuscript, which defines the rational series resolution:

$$\sum_{k=0}^\infty \frac{P(k)\, z^k}{(k+a_1)^{m_1+1} (k+a_2)^{m_2+1} \cdots (k+a_n)^{m_n+1}} = \sum_{i=1}^{n} \frac{1}{\prod_{\substack{j=1 \\ j \ne i}}^{n} (a_j - a_i)^{m_j + 1}} \sum_{k=0}^{m_i} \frac{(-1)^kF^{(k)}(z, a_i)}{k!} \, \Omega_{i,\, m_i - k}$$

Under specific adaptations, this theoretical development directly translates to the physical problem of nuclear decay chains. The general solution of the Bateman equations with repeated decay constants shares the exact same algebraic-combinatorial structure, as shown in **Equation (NewC_01)**:

$$\frac{X_N(t)\lambda_N}{X_1(0)\prod_{k=1}^{n}\lambda_k^{\mu_k+1}} = \sum_{i=1}^{n} \frac{1}{\prod_{\substack{j=1 \\ j \neq i}}^{n} (\lambda_j - \lambda_i)^{\mu_j}} \sum_{l=0}^{\mu_i} \frac{(-1)^l G^{(l)}(\lambda_i,t)}{l!}\,\chi_{i,\mu_i-l}$$

This structural corollary implies that the algorithm developed for evaluating the rational series (Eq. 60) is directly adapted to compute the general Bateman solution.

---

# Computational Implementation

The numerical evaluation of these analytical solutions requires high-precision computations. The scripts are implemented in **Python 3**, utilizing:
* `mpmath` for arbitrary-precision floating-point arithmetic,
* `sympy` for symbolic mathematics,
  
# Symbolic Verification and PSLQ Analysis

To ensure the absolute theoretical rigidity and exactness of the proposed analytical solutions, this repository also includes scripts developed in **Wolfram Mathematica**. 

These scripts complement the numerical framework by providing:
* **Symbolic Computation:** Exact algebraic verification of the rational series resolutions, the polynomial coefficients, and the general structure of the Bateman decay chains.
* **PSLQ Algorithm Implementation:** Application of the high-precision PSLQ integer relation algorithm to rigorously validate the analytical structure, detect precise integer relations among the coefficients, and confirm the exactness of the developed framework.
* 
## Acknowledgments

The authors Carlos Antonio Cruz López and Gilberto Espinosa Paredes gratefully acknowledge the financial support received from the Secretaría de Ciencia, Humanidades, Tecnología e Innovación (SECIHTI, formerly known as CONAHCYT), through the program *Estancias Posdoctorales por México, 2022*, under the project entitled: *Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia*, as well as the financial support from the Basic Science and Frontier Project 2023-2024, with the reference CBF-2023-2024-2023, also belonging to SECIHTI-CONAHCYT. 

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. You are free to use, modify, and distribute this software for academic and commercial purposes, provided that appropriate credit is given to the original authors.
