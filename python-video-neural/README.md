# Neuronal responses to natural scenes with modified spatial correlations.

Refer paper for information about experiment and setup.

## About data:
MATLAB datafile AmpMov.mat contains the following fields.

Experiments are done various days, the data corresponding to each day are in each folder.

Subscript _nat corresponds to responses to video stimuli for original natural scenes video. similarly subscripts _K0, _K_1, _K1_5, K_2 denote responses to manipulated natural scenes video stimuli.

Duration of movies = 4s at 20Hz with a blank screen for 2 seconds.
Field                               | Description 
---                                 | -------
1. Sorted.SpikeRate                 | blah blah
2. Blank                            | Dimension 47 x 16800
2. NumNeurons                       | Number of neurons sampled
3. NumMovies                        | Number of movies used as stimulus

4. M_nat                            | Dimension 4 x 47 x 1200 
5. MT_nat                           | Dimension 4 x 47 x 200 x 6
6. MTA_nat                          | Dimension 4 x 47 x 200
7. MTNA_nat                         | Dimension 4 x 47 x 1200