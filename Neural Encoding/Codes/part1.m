%%
clc
clear all
close all
T = load('UnitsData.mat');
%%
M = T.Unit;
%%
numOfTrials = 192;
numOfUnits = 481;
numOfBins = 32;
meanOfFreq1 = zeros(numOfUnits, 32);
meanOfFreq2 = zeros(numOfUnits, 32);
meanOfFreq3 = zeros(numOfUnits, 32);
meanOfFreq4 = zeros(numOfUnits, 32);
meanOfFreq5 = zeros(numOfUnits, 32);
meanOfFreq6 = zeros(numOfUnits, 32);
for k = 1 : numOfUnits
    for i = 1 : numOfTrials
        temp = M(k).Trls{i, 1};
        s = zeros(1, numOfBins);
        for j = 1 : numOfBins
            s(j) =  length(temp(temp>= (j-1)*3.2/numOfBins-1.2 & temp<=(j)*3.2/numOfBins-1.2));
        end
        if(ismember(i,M(k).Cnd(1).TrialIdx))
            meanOfFreq1(k, :) = meanOfFreq1(k,:) + s./length(M(k).Cnd(1).TrialIdx);
        elseif(ismember(i,M(k).Cnd(2).TrialIdx))
            meanOfFreq2(k, :) = meanOfFreq2(k,:) + s./length(M(k).Cnd(2).TrialIdx);
        elseif(ismember(i,M(k).Cnd(3).TrialIdx))
            meanOfFreq3(k, :) = meanOfFreq3(k,:) + s./length(M(k).Cnd(3).TrialIdx);
        elseif(ismember(i,M(k).Cnd(4).TrialIdx))
            meanOfFreq4(k, :) = meanOfFreq4(k,:) + s./length(M(k).Cnd(4).TrialIdx);
        elseif(ismember(i,M(k).Cnd(5).TrialIdx))
            meanOfFreq5(k, :) = meanOfFreq5(k,:) + s./length(M(k).Cnd(5).TrialIdx);
        elseif(ismember(i,M(k).Cnd(6).TrialIdx))
            meanOfFreq6(k, :) = meanOfFreq6(k,:) + s./length(M(k).Cnd(6).TrialIdx);
        end
    end
end
%% Average PSTH
close all
figure;
plotPSTH(mean(meanOfFreq1)); hold on;
plotPSTH(mean(meanOfFreq2)); hold on;
plotPSTH(mean(meanOfFreq3)); hold on;
plotPSTH(mean(meanOfFreq4)); hold on;
plotPSTH(mean(meanOfFreq5)); hold on;
plotPSTH(mean(meanOfFreq6)); 
legend('CND 1', 'CND 2', 'CND 3', 'CND 4', 'CND 5', 'CND 6', 'Location','southwest')
xlabel("time( aligned to Cue Onset) sec")
ylabel("firing rate")
title("average PSTH for each condition of the task over all units")
hold off;
%% pSTH for different neurons
close all
for i = 241 : 241
    plotPSTH1Unit(meanOfFreq1, meanOfFreq2, meanOfFreq3, meanOfFreq4, meanOfFreq5, meanOfFreq6, i)
%     saveas(gcf, sprintf("%d.png", i))
end
%% regression
Stats = zeros(481, 4);
regCoeff = zeros(481, 3);
ys = zeros(481, 6);
X = [3 1 1; 3 -1 1; 6 1 1; 6 -1 1; 9 1 1; 9 -1 1];
for i = 1 : 481
ys(i, :) = [mean(10*meanOfFreq1(i,16:26))
    mean(10*meanOfFreq2(i,16:26))
    mean(10*meanOfFreq3(i,16:26))
    mean(10*meanOfFreq4(i,16:26))
    mean(10*meanOfFreq5(i,16:26))
    mean(10*meanOfFreq6(i,16:26))];
[b,~,~,~,stats] = regress(ys(i, :)',X);
Stats(i, :) = stats;
regCoeff(i, :) = b;
end
%%
CV = zeros(481, 1);
for i = 1 : 481
    CV(i) = sqrt(var(ys(i, :)))/abs(mean(ys(i, :)));
end

%%
[cvsort, idx] = sort(CV);
CVUnits = idx(332:end);
chosenUnits = [CVUnits Stats(CVUnits, 3)];
sortedChosenUnits = sortrows(chosenUnits,2);
%%
impIdx = sortedChosenUnits([1:5:26 26:25:126] , 1);
%%
close all
for i = 1:11
    plotRegression(X, ys(impIdx(i),:), regCoeff(impIdx(i),:), impIdx(i), Stats(impIdx(i), [1 3]), CV(impIdx(i)))
    fname = 'E:\UNI\advanced neuro\HW2\fig';
saveas(gca, fullfile(fname, sprintf("%d.png", i+19)))
end
%% pca
close all
y = [mean(meanOfFreq1(:, 16:26)')' mean(meanOfFreq2(:,16:26)')' mean(meanOfFreq3(:,16:26)')' mean(meanOfFreq4(:,16:26)')' mean(meanOfFreq5(:,16:26)')' mean(meanOfFreq6(:,16:26)')'];
[coeff,score,latent,tsquared,explained,mu] = pca(y)
% scatter(score(:, 4), score(:, 2))
%% main pca
close all
y = meanOfFreq1(:, 16:26);
[coeff,score,latent,tsquared,explained,mu] = pca(y);

y1 = score(:, 3);
y2 = score(:, 2);
t = 1:200;
plot(y1(t), y2(t))




%% functions

function plotPSTH(meanOfFreq)
    x = [1:1000]./312.5 - 1.2;
    y = zeros(1,1000);
    siz = size(meanOfFreq(:, :));
    yTemp = meanOfFreq;
    for j = 1 : siz(2)
        y(floor(1000/32*(j-1))+1 : floor(1000/32*j)) = yTemp(j);
    end
    plot(x, 10*y)
end

function plotPSTH1Unit(meanOfFreq1, meanOfFreq2, meanOfFreq3, meanOfFreq4, meanOfFreq5, meanOfFreq6, k)
figure;
plotPSTH(meanOfFreq1(k, :)); hold on;
plotPSTH(meanOfFreq2(k, :)); hold on;
plotPSTH(meanOfFreq3(k, :)); hold on;
plotPSTH(meanOfFreq4(k, :)); hold on;
plotPSTH(meanOfFreq5(k, :)); hold on;
plotPSTH(meanOfFreq6(k, :)); 
legend('CND 1', 'CND 2', 'CND 3', 'CND 4', 'CND 5', 'CND 6', 'Location','southwest')
xlabel("time( aligned to Cue Onset) sec")
ylabel("firing rate")
title(sprintf("PSTH for each condition of the task over corresponding trials of Unit %d", k))
hold off;
end

function plotRegression(X, y, b, k, m, n)
f = figure;
f.Position = [480   258 680 510];
x1 = X(:, 1);
x2 = X(:, 2);
scatter3(x1,x2,y,'filled')
hold on
x1fit = 0:0.05:10;
x2fit = -2:0.05:2;
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(3) + b(1)*X1FIT + b(2)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Expected Reward')
ylabel('Cue Location')
zlabel('Firing rate')
title(["regressed neural responses against cue location and expected reward", sprintf("for unit %d with P-Value = %f and R-squared = %f", k, m(2), m(1)), sprintf("CV = %f", n)])
view(70,20)
hold off

end



