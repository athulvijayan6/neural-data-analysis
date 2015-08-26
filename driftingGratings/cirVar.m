% @Author: athul
% @Date:   2015-08-26 02:14:07
% @Last Modified by:   athul
% @Last Modified time: 2015-08-26 02:30:43
%% cirVar: function description
function [cv] = cirVar(spikeRate)
    cv = (transpose(spikeRate(:, 1)) * exp(2i*spikeRate(:, 2)))/sum(spikeRate(:, 1));
end
    
