---
layout: post
title:  "Test Signals (ssqueezepy)"
date:   2021-01-30 00:01:10 +0400
author: OverLordGoldDragon
categories: jekyll update
---

### Test signals

To gauge tradeoffs between time and frequency localization, we design signals varying in each; no single
wavelet or transform can excel at accurately mapping all.

 - `#`: reversed signal is added to itself (`x += x[::-1]`)
 - `amin`: amplitude modulator minimum (default `amax=1`)
 - `lchirp, echirp, hchirp`: linear, exponential, hyperbolic chirp

Code at [ssqueezepy](https://github.com/OverLordGoldDragon/ssqueezepy/blob/master/examples/test_transforms.py).

<img src="https://user-images.githubusercontent.com/16495490/106352545-84d91800-62fd-11eb-900e-4703178a5246.png">

### Wavelets vs test signals; CWT

Compare Generalized Morse Wavelets with `beta=5` (high time localization) vs `beta=22` (high frequency localization), `gamma=3` for both (optimal join localization).

<img src="https://user-images.githubusercontent.com/16495490/106352592-d4b7df00-62fd-11eb-8ccc-21d456541751.png">


### Signal general forms derivations

Python pseudocode

```python
"""
---- **lchirp** --------------------------------------------------------------
https://www.desmos.com/calculator/cqe5oqlnyt
Design LFM to range from `(tmin, fmin)` to `(tmax, fmax)`

Choose form: phi(t) = A*t^2 + B*t + K --> f(t) = a*t + b

>> (tmin, fmin)
f(tmin) = a*tmin + b = fmin
-> b = fmin - a*tmin

>> (tmax, fmax)
f(tmax) = a*tmax + b = fmax
-> b = fmax - a*tmax

-> fmax - a*tmax = fmin - a*tmin
   a*(tmin - tmax) = fmin - fmax
   a = (fmin - fmax) / (tmin - tmax) ####

-> b = fmin - a*tmin
     = fmin - tmin * (fmin - fmax)/(tmin - tmax)
     = N / D
> D = (tmin - tmax)
> N = fmin*(tmin - tmax) - tmin*(fmin - fmax)
    = fmin*tmin - fmin*tmax - tmin*fmin + tmin*fmax
    = fmin*tmin - fmin*tmax - fmin*tmin + fmax*tmin
    = fmax*tmin - fmin*tmax
> N / D = (fmax*tmin - fmin*tmax) / (tmin - tmax)
        = (fmin*tmax - fmax*tmin) / (tmax - tmin)
        = b ####

>> phi(t)
phi(t) = int_{tmin}^{t} f(t) dt
       = [(a/2)*t^2 + b*t + C]_{tmin}^{t}
       = [(a/2)*t^2 + b*t] - [(a/2)*tmin^2 + b*tmin]
       = (a/2)*(t^2 - tmin^2) + b*(t - tmin) ####

>>> f(t)   = a*t + b
>>> phi(t) = (a/2)*(t^2 - tmin^2) + b*(t - tmin)
>>> a = (fmin - fmax) / (tmin - tmax)
    b = (fmin*tmax - fmax*tmin) / (tmax - tmin)


^ simple example for general procedure. Simpler derivation, 
  use point-slope form:
    y - y0 = m*(x - x0)
    y = mx + y0 - m*x0, m = (y1 - y0) / (x1 - x0)

f(t) = m*(t - tmin) + fmin, m = (fmax - fmin) / (tmax - tmin)
     = m*t + (fmin - m*tmin)
phi(t) = [(m/2)*t^2 + (fmin - m*tmin)*t]_{tmin}^{t}
       = [(m/2)*t^2 + (fmin - m*tmin)*t] -
         [(m/2)*tmin^2 + (fmin - m*tmin)*tmin]
       = (m/2)*(t^2 - tmin^2) + (fmin - m*tmin)*(t - tmin)


---- **echirp** --------------------------------------------------------------
https://www.desmos.com/calculator/jaile5fabs
Design EFM to range from `(tmin, fmin)` to `(tmax, fmax)`

Choose form: phi(t) = A*exp(t) + B*T + K --> f(t) = a*exp(t) + b
This keeps f'(t) independent of `tmin, fmin, tmax, fmax`.

>> (tmin, fmin)
f(tmin) = a*exp(tmin) + b = fmin
-> b = fmin - a*exp(tmin)

>> (tmax, fmax)
f(tmax) = a*exp(tmax) + b = fmax
-> b = fmax - a*exp(tmax)

-> fmax - a*exp(tmax) = fmin - a*exp(tmin)
   a*(exp(tmin) - exp(tmax)) = fmin - fmax
   a = (fmax - fmin)/(exp(tmax) - exp(tmin)) ####

-> b = fmin - a*exp(tmin)
     = fmin - exp(tmin) * (fmax - fmin)/(exp(tmax) - exp(tmin))
     = N / D
> D = (exp(tmax) - exp(tmin))
> N = fmin*(exp(tmax) - exp(tmin)) - exp(tmin)*(fmax - fmin)
    = fmin*exp(tmax) - fmin*exp(tmin) - fmax*exp(tmin) + fmin*exp(tmin)
    = fmin*exp(tmax) - fmax*exp(tmin)
> N / D = (fmin*exp(tmax) - fmax*exp(tmin)) / (exp(tmax) - exp(tmin))
        = b ####

>> phi(t)
phi(t) = int_{tmin}^{t} f(t) dt
       = [a*exp(t) + b*t + C]_{tmin}^{t}
       = [a*exp(t) + b*t] - [a*exp(tmin) - b*tmin]
       = a*(exp(t) - exp(tmin)) + b*(t - tmin) ####

>>> f(t)   = a*exp(t) + b
>>> phi(t) = a*(exp(t) - exp(tmin)) + b*(t - tmin)
>>> a = (fmax - fmin)/(exp(tmax) - exp(tmin))
    b = (fmin*exp(tmax) - fmax*exp(tmin)) / (exp(tmax) - exp(tmin))


---- **hchirp** --------------------------------------------------------------
https://www.desmos.com/calculator/v1iu9ydjrq
Design HFM to range from `(tmin, fmin)` to `(tmax, fmax)`

Choose form: phi(t) = a / (b - t) --> f(t) = A / (B - t)^2

>>> f(t)   = A / (B - t)^2
>>> phi(t) = A * (1/(B - t) + 1/(tmin - B))
>>> a, b, c, d = fmin, fmax, tmin, tmax
    A = AN / AD, B = BN / BD,
    AN = 2*sqrt(a^3*b^3*(c - d)^4) + a^2*b*(c - d)^2 + a*b^2*(c - d)^2
    AD = (a - b)^2
    BN = sqrt(a^3*b^3*(c-d)^4) + a^2*b*c*(c-d) + a*b^2*d*(d - c)
    BD = a*b*(a - b)*(c - d)

Derivation: tedious, feed to wolf instead:
https://www.wolframalpha.com/input/
?i=solve+x+%3D+a*%28y+-+c%29%5E2%2C+x+%3D+b*%28y+-+d%29%5E2
"""
```
