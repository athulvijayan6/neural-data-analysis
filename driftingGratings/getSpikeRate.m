% @Author: Athul Vijayan
% @Date:   2015-08-22 12:51:54
% @Last Modified by:   athul
% @Last Modified time: 2015-08-26 02:20:46

clear('all');

datTargets = {'dataset/Mouse-A/', 'dataset/Mouse-B/', 'dataset/Mouse-C/', 'dataset/Mouse-D/', 'dataset/Mouse-E/'};

% for each mouse
load('dataset/Mouse-A/Data.mat');
rawData = Data.rawF;
smoothData = Data.dFF;
spikeData = Data.Spks;
stimuliSeq = Data.StimSeq;
numClasses = 16;
for i=1:size(rawData, 1)
    for j=1:size(stimuliSeq, 1)
        cellData{i}(j, :) = [smoothData(i, 120*(j-1)+1:120*j), stimuliSeq(j)];
    end
    % sort the array
    [~, I] = sort(cellData{i}(:, end));
    cellData{i} = cellData{i}(I, :);

    spikeRate{i} = calculateSpikeRate(cellData{i});
    s{i} = cirVar(spikeRate{i});
end

% ====================== Plots =======================
% hold on;
% Plots each time-series data of length 120 
% for j=1:size(stimuliSeq, 1)
%     plot(1:120, cellData{1}(j, 1:120));
% end

% plot(1:19200, reshape(cellData{1}(:, 1:120), 1, prod(size(cellData{1}(:, 1:120)))));


% % plots The responses with increasing angle including the gray region
% plot(1:1920, reshape(cellData{1}(1:10:160, 1:120), 1, prod(size(cellData{1}(1:10:160, 1:120)))))

% % plots The responses with increasing angle without the gray region
% plot(1:1920, reshape(cellData{1}(1:10:160, 1:120), 1, prod(size(cellData{1}(1:10:160, 1:120)))))

% ********************** Plots ************************