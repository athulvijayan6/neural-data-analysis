### Module for orientation and directional selectivity of neurons

## Methods
---------
1.`calculateSpikeRate(conc, algo='average')`: the function computes the neuron response of a stimulus using the time-series. The function returns a real value which can be thought of as the spike rate during the stimulus period.
2.`cirVar(spikeRate)`: Function calculates the L_ori from spike rate response of neuron to each angle and direction of stimuli
3.`dirCirVar(spikeRate)`: Function calculates the L_dir from spike rate response of neuron to each angle and direction of stimuli
