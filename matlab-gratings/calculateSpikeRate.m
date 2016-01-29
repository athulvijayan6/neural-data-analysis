% @Author: athul
% @Date:   2015-08-26 02:03:56
% @Last Modified by:   Athul
% @Last Modified time: 2015-09-01 21:31:46

% ======================= Functions ==================
%% getSpikeRate: function description
function [spikeRate] = calculateSpikeRate(conc, varargin)
    if nargin > 1
        algo = varargin{1};
    else
        algo = 'average';
    end
    spikeRate = zeros(size(conc, 1), 2);
    for i=1:size(conc, 1)
        if algo == 'average'
            % This algorithm just computes the average Ca concentration
            % from the time of epoch to end of stimuli.The input expected is smoothed Ca data
            spikeRate(i, 1) = mean(conc(i, 81:end));
            spikeRate(i, 2) = conc(i, end);
        end
    end
end