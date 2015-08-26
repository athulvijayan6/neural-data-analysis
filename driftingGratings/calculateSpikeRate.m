% @Author: athul
% @Date:   2015-08-26 02:03:56
% @Last Modified by:   athul
% @Last Modified time: 2015-08-26 02:08:15

% ======================= Functions ==================
%% getSpikeRate: function description
function [spikeRate] = calculateSpikeRate(conc, varargin)
    if nargin > 1
        algo = varargin{1};
    else
        algo = 'threshold';
    end
    spikeRate = zeros(size(conc, 1), 2);
    for i=1:size(conc, 1)
        if algo == 'threshold'
            spikeRate(i, 1) = max(conc(i, 1:end-1));
            spikeRate(i, 2) = conc(i, end);
        end
    end
end