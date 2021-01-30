---
layout: post
title:  "Generalized Morse Wavelets"
date:   2021-01-19 00:01:10 +0400
author: OverLordGoldDragon
categories: jekyll update
---

GMWs are defined in the frequency domain, as:

$$
\Psi_{k; \beta, \gamma}(\omega^+) 
= A_{k; \beta, \gamma} \omega^\beta e^{-\omega^\gamma} L_k^c(2 \omega^\gamma)
$$

where $L_k^c$ denotes the generalized Laguerre polynomial

$$
L_k^c (x)
 = \sum_{m=0}^k (-1)^m 
 \frac{ \Gamma(k + c + 1) }{ \Gamma(c + m + 1)\Gamma(k - m + 1) }
 \frac{x^m}{m!}
$$

where $c = r - 1$ and $r = (2\beta + 1) / \gamma$, and the amplitude term $A$ depends on normalization, L1 (bandpass) or L2 (energy): 

$$
\begin{align}
A_{k; \beta, \gamma}^{\text{L1}} 
 &= 2 \left( \frac{e \gamma}{\beta} \right)^{\beta / \gamma} \\
A_{k; \beta, \gamma}^{\text{L2}}
 &= \sqrt{ 2\pi \gamma 2^r \frac{\Gamma(k + 1)}{\Gamma(k + r)}
           \left( \frac{\beta}{\gamma} \right)^{1 / \gamma}
          } \\
\end{align}
$$

$\omega^+$ replaces $U(\omega)$ (unit-step function), for simplicity, in denoting $\Psi(\omega < 0) = 0$ (analyticity).

<hr>

In single expression:

$$
\begin{align}
\Psi_{k; \beta, \gamma}^{\text{L1}}(\omega^+) 
 &= 2 \left( \frac{e \gamma}{\beta} \right)^{\beta / \gamma}
    \omega^\beta e^{-\omega^\gamma} 
    \sum_{m=0}^k (-1)^m 
    \frac{ \Gamma(k + r) }{ \Gamma(r + m)\Gamma(k - m + 1) }
    \frac{(2\omega^\gamma)^m}{m!} \\
\Psi_{k; \beta, \gamma}^{\text{L2}}(\omega^+) 
 &= \sqrt{ 2\pi \gamma 2^r \frac{\Gamma(k + 1)}{\Gamma(k + r)}
          \left( \frac{\beta}{\gamma} \right)^{1 / \gamma}
          }
 \omega^\beta e^{-\omega^\gamma} 
 \sum_{m=0}^k (-1)^m 
 \frac{ \Gamma(k + r) }{ \Gamma(r + m)\Gamma(k - m + 1) }
 \frac{(2\omega^\gamma)^m}{m!} \\
\end{align}
$$

For the base wavelet, where $k=0$ and $L_0^c = 1$, these simplify to:

$$
\begin{align}
\Psi_{\beta, \gamma}^{\text{L1}} (\omega^+)
 &= 2 \left( \frac{e \gamma}{\beta} \right)^{\beta / \gamma}
    \omega^\beta e^{-\omega^\gamma} \\
\Psi_{\beta, \gamma}^{\text{L2}} (\omega^+)
 &= \sqrt{ 2\pi \gamma 2^r \frac{1}{\Gamma(r)} 
           \left( \frac{\beta}{\gamma} \right)^{1 / \gamma} 
          }
 \omega^\beta e^{-\omega^\gamma} \\
\end{align}
$$

and the (peak) center frequency for both is $\omega_{\beta, \gamma} = (\beta / \gamma)^{1/\gamma}$.

<hr>

### Notes

**Scale**: above expressions omit _scale_; when accounted, L1 is computed as-is, simply
inputting $\omega := s \cdot \omega$, but L2 additionally multiplies by $\sqrt{s}$:

$$
\begin{align}
\Psi_{s, k; \beta, \gamma}^{\text{L1}} (\omega)
 &= \Psi_{k; \beta, \gamma}^{\text{L1}} (s \omega) \\
\Psi_{s, k; \beta, \gamma}^{\text{L2}} (\omega)
 &= \sqrt{s} \Psi_{k; \beta, \gamma}^{\text{L2}} (s \omega) \\
\end{align}
$$

**Analytic vs anti-analytic**: above expressions are for the _analytic_ GMW, $\Psi^+$; 
there's also the _anti-analytic_ $\Psi^-$, non-zero for _negative_ frequencies, which 
computes its conjugate in time domain: 
$\Psi^- (\omega) = \Psi^+ (-\omega) \Leftrightarrow \psi^-(t) = \psi^{+*}(t)$. 
Combining the two thus yields a real wavelet.

**Interactive plots**: [GMW L1](https://www.desmos.com/calculator/s3mam3qhp3) -- 
[GMW L2](https://www.desmos.com/calculator/lttewqu89w)

<hr>

### Derivations

Python pseudocode deriving above per original (MATLAB, JLAB) formulation (w/ name changes), 
following 
[`.m` code](https://github.com/jonathanlilly/jLab/blob/master/jWavelet/morsewave.m), 
for the case `k=0` where `L=1`:

```python
# note: `f` and `f0` are *radian* quantities, not linear cyclic
# `**` and `^` used interchaneably
f0 = exp( (ln(beta) - ln(gamma)) / gamma )
   = exp( ln(beta^{1/gamma}) - ln(gamma^{1/gamma}) )
   = beta^{1/gamma} / gamma^{1/gamma}
   = (beta / gamma)^{1/gamma}

## norm='bandpass' ############################################
psizero = 2 * exp(- beta * ln(f0) + f0**gamma
                  + beta * ln(w)  - w**gamma)
        = 2 * H * G 

H = exp( beta * (ln(w) - ln(f0)) )
  = exp( beta * (ln(w) - ln((beta/gamma)^{1/gamma})) )
  = exp( ln(w^beta) - ln((beta/gamma)^{beta/gamma}) )
  = exp( ln(w^beta) + ln((gamma/beta)^{beta/gamma}) )
  = w^beta (gamma/beta)^{beta/gamma}

G = exp( fo^gamma - w^gamma )
    exp( (beta/gamma)^(gamma/gamma) - w^gamma )
    exp( beta/gamma - w^gamma )

Psi = psizero  # coeff == 1, L == 1
    = 2 * H * G
    = 2 * w^beta * (gamma/beta)^{beta/gamma} * e^{ beta/gamma - w^gamma }
    = 2 * w^beta * (gamma*e/beta)^{beta/gamma} * e^{-w^gamma}
    = 2 * (gamma*e/beta)^{beta/gamma} * w^beta * e^{-w^gamma}
       
## norm='energy' ##############################################
psizero = exp( beta * ln(w) - w**gamma )
        = exp( ln(w^beta) - w^gamma )
        = w^beta * e^{-w^gamma}  

A = sqrt( 2*pi * gamma * (2**r) * exp(gammaln(k+1) - gammaln(k + r)) )
  = sqrt( 2*pi * gamma * (2**r) * exp( -ln(gamma(r)) ) )
  = sqrt( 2*pi * gamma * (2**r) * exp( ln(1/gamma(r)) ) )
  = sqrt( 2*pi * gamma * (2**r) / gamma(r) )

coeff = sqrt(f0/f) * A
      = sqrt( 2*pi * (f0/f) * gamma * (2**r) / gamma(r) )

Psi = coeff * psizero * L
    = sqrt( 2*pi * gamma * (2^r) / gamma(r) ) * sqrt(f0/f)
      * w^beta * e^{-w^gamma}
    = sqrt( 2*pi * gamma * (2^r) / gamma(r) * (beta / gamma)^{1/gamma} / f)
      * w^beta * e^{-w^gamma}
```

Note that `sqrt(f0/f)` is absent in Eq 10 of ref 1 (L2 norm case). Its inclusion 
normalizes the time-domain wavelet's energy to unity, and the complete expression
is given in CWT examples as `sqrt(scale) * psih`.

<hr>

### References

 1. ["Generalized Morse Wavelets".](https://spiral.imperial.ac.uk/bitstream/10044/1/1150/1/OlhedeWaldenGenMorse.pdf)
 Sofia C. Olhede, Andrew T. Walden. November, 2002. <br>_IEEE Transactions On Signal Processing_.
 2. ["Higher-Order Properties of Analytic Wavelets".](https://sci-hub.st/10.1109/TSP.2008.2007607)
 Johnathan M. Lilly, Sofia C. Olhede. January, 2009. _IEEE Transactions On Signal Processing_.
