% @Author: Athul Vijayan
% @Date:   2015-08-22 12:51:54
% @Last Modified by:   Athul Vijayan
% @Last Modified time: 2015-08-29 14:22:51

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
colors = hsv(10);
neuronId = 20;

% -------------------------1------------------------------
% % 1. Plots each time-series data of length 120 with 
% % each trial of an experiment with seperate color
% figure;
% hold on;
% experimentNo = 2;
% for j=1:10
%     h = plot(1:120, cellData{1}(10*(experimentNo-1) + j, 1:120));
%     set(h,'Color', colors(j, :));
%     set(h,'LineWidth',2);
% end
% title('Plots of time-series data with each trial in diffenrent color');
% xlabel('time');
% ylabel('smoothed amplitude');
% --------------------------------------------------------

% --------------------------2-----------------------------
% 2. plots The responses with increasing angle without the gray region
% plot(1:1920, reshape(cellData{1}(1:10:160, 1:120), 1, prod(size(cellData{1}(1:10:160, 1:120)))));
% ---------------------------------------------------------

% --------------------------3-----------------------------
% 3. Plots the gray area alone. To show that the background data are almost uniform
% --------------------------------------------------------


% --------------------------4-----------------------------
% 4. Polar plot of Directional selectivity
% Take a neuron, there will be 10 plots corresponding to each trial of the experiment
% plot each of the 10 plots in a figure to verify the similarity in response
figure;
for trial=1:10
    % Normalize all trial response from 0 to 1
    normData = abs(spikeRate{neuronId}(trial:10:end, 1))/max(abs(spikeRate{neuronId}(trial:10:end, 1)));
    normData(end+1) = normData(1);
    theta = spikeRate{neuronId}(trial:10:end, 2);
    theta(end+1) = theta(1);
    h = polar(theta, normData);
    set(h,'Color',colors(trial, :));
    set(h,'LineStyle','--');
    set(h,'Marker','o');
    set(h,'LineWidth',2);
    hold on;
end
h = compass(dircirvar{neuronId});
set(h,'LineWidth',5);
set(h,'Color', 'r');
title('Polar plot of Directional selectivity');
% --------------------------------------------------------

% --------------------------5-----------------------------
% 4. Polar plot of Orientation selectivity
% Take a neuron, there will be 10 plots corresponding to each trial of the experiment
% plot each of the 10 plots in a figure to verify the similarity in response

% There are observations for 16 angles, but there are only 8 orientations.
% We can treat them as two different observations. More data !

figure;
for trial=1:10
    % Normalize all trial response from 0 to 1
    normData = abs(spikeRate{neuronId}(trial:10:end, 1))/max(abs(spikeRate{neuronId}(trial:10:end, 1)));
    normData(end+1) = normData(1);
    theta = 2*spikeRate{neuronId}(trial:10:end, 2);
    theta(end+1) = theta(1);
    h = polar(theta, normData);
    set(h,'Color',colors(trial, :));
    set(h,'LineStyle','--');
    set(h,'Marker','o');
    set(h,'LineWidth',2);
    hold on;
end
h = compass(cirvar{neuronId});
set(h,'LineWidth',5);
set(h,'Color', 'r');
title('Polar plot of Orientation selectivity');
% --------------------------2-----------------------------

% ********************** Plots ************************