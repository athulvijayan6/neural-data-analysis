% @Author: athul
% @Date:   2015-08-26 02:15:08
% @Last Modified by:   Athul
% @Last Modified time: 2015-08-27 20:55:40

%% dirCirVar: Computes the directional circular variance.
function [dcv] = dirCirVar(spikeRate)
    dcv = (transpose(spikeRate(:, 1)) * exp(i*spikeRate(:, 2)))/sum(spikeRate(:, 1));
end
