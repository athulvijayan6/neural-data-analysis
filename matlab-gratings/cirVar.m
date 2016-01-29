% @Author: athul
% @Date:   2015-08-26 02:14:07
% @Last Modified by:   Athul
% @Last Modified time: 2015-09-01 21:31:57
%% cirVar: function description
function [cv, points] = cirVar(spikeRate)
    cv = (transpose(spikeRate(:, 1)) * exp(2i*spikeRate(:, 2)))/sum(spikeRate(:, 1));
end
    
