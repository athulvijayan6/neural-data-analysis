% @Author: Athul Vijayan
% @Date:   2015-08-22 12:51:54
% @Last Modified by:   Athul
% @Last Modified time: 2015-08-28 13:14:36

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
        cellData{i}(j, :) = [smoothData(i, 120*(j-1)+1:120*j), degtorad(stimuliSeq(j))];
    end
    % sort the array
    [~, I] = sort(cellData{i}(:, end));
    cellData{i} = cellData{i}(I, :);

    % Calculate response of each neuron as a real number
    spikeRate{i} = calculateSpikeRate(cellData{i});
    cirvar{i} = cirVar(spikeRate{i});
    dircirvar{i} = dirCirVar(spikeRate{i});
end

% ====================== Plots =======================
% 1. Plots each time-series data of length 120 
    % figure;
    % hold on;
    % for j=1:size(stimuliSeq, 1)
    %     plot(1:120, cellData{1}(j, 1:120));
    % end

% 2. plots The responses with increasing angle without the gray region
    % plot(1:1920, reshape(cellData{1}(1:10:160, 1:120), 1, prod(size(cellData{1}(1:10:160, 1:120)))))

% 3. Plots the gray area alone. To show that the background data are almost uniform

% 4. Polar plot of orientation selectivity
figure
% Take a neuron, there will be 10 plots corresponding to each trial of the experiment
% plot each of the 10 plots in a figure to verify the similarity in response
neuronId = 20;
colors = ['']
for 
polar(spikeRate{1}(1, 1), spikeRate{1}(1, 2), '*r');

% ********************** Plots ************************