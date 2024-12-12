# mRNoise
## Radio noise model for Matlab based on ITU-R P.372 Recommendation

This package provides worldwide radio noise calculations. The core of the package consists of several interpolant objects. This interpolators are saved in the Noise.mat file. All functions on the lowest level load interpolators from the Noise.mat file. But it is possible to use them directrly.

There are two types of functions. 

First type: `DParam1`, `NoiseFam`, `NoiseTotal` and `NoisePwr` calculate spatial noise distribution.

Second Type: `DParam`, `NoiseFamPoint`, `NoiseTotalPoint` and `NoisePwrPoint` calculate diurnal variations.

Function `NoiseMG` is universal because it is not tied to time or coordinates.

The main functions are `NoiseTotal` and `NoisePwr`. `NoiseTotal` returns structure with median total noise figure, median values of noise figure components and appropriate upper and lower deciles. `NoisePwr` returns noise power and E-field. During calculations `NoiseTotal` and `NoisePwr` call other functions and interpolators. Also `NoisePwr` calls `NoiseTotal`.

All inputs and outputs described within approprite functions files. Script `Noise_examples.m` contains usage examples. Before running change `fpath` variable to the appropriate path (for functions) or `Noise.mat` full name (for interpolators).
