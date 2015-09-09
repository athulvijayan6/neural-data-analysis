# Calcium signals (GCaMP6f in awake mice) for visual stimuli of Sinusoidal drifting gratings

## About the structure of data
-----

### Notes about Data.mat

Data.mat contains 4 entries:
1. Data.rawF =  raw fluorescence values. Matrix size = number of cells x number of frames.
2. Data.dFF    = fluorescence normalized to baseline (dFF =  (F-F0)/F0, where F0 is the baseline fluorescence computed using a sliding window of 400 frames). Same size as above.
3. Data.Spks  = inferred spike rate using the Vogelstein deconvolution algorithm. Same size as above.
4. Data.StimSeq = contains sequence of directions presented during that experiment. Vector size = 160 x 1.
*****
### Notes about Ori.mat

Ori.mat contains 20 entries. The most pertinent entries are:
1. Ori.OSI = orientation selectivity indices of the neurons
2. Ori.OrFit = double-wrapped Gaussian fits 
3. Ori.OrFitQuality  = goodness of Gaussian fits. (Higher the percentage value, the better the fit)
4. Ori.Width = tuning width in degrees
5. Ori.PrefOri = preferred orientation 
6. Ori.SpkResponse = Contains neural responses for each cells sorted according to the different directions. Size: 1xNumber of cell Cell array. Each cell entry contains a 1x Number of Direction Cell array, which contains a Number of Frames x Trials matrix.

*****