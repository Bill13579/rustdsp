# rustdsp
The ultimate audio filters library written for Rustlang.

## Generating various filters using filter transformations

Various filter types can be achieved through the use of these, for example, a high pass can be achieved by the use of `InvertCutoff`

## Overview on the Fractional Coefficients

Fractional coefficients were generated for the entire set of slopes from 6-120 dB/oct using Dual Annealing on the set of parameters for a set of biquads ( biquads), with stability ensured through restricting the position of the poles in the biquads to the left side of the complex S-plane, to the negative real side.
More detailed information will be posted as development continues.

## Stability

This library is still in **heavy development**, all the filters are working, but things might change at any time.
