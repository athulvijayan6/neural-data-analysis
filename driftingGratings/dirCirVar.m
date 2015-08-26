% @Author: athul
% @Date:   2015-08-26 02:15:08
% @Last Modified by:   athul
% @Last Modified time: 2015-08-26 02:31:50

%% dirCirVar: Computes the directional circular variance.
function [dcv] = dirCirVar(spikeRate)
    cv = (transpose(spikeRate(:, 1)) * exp(i*spikeRate(:, 2)))/sum(spikeRate(:, 1));
end
