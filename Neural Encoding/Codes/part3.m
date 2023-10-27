clc
clear all
close all
T = load('UnitsData.mat');
%%
M = T.Unit;
%% shuffle
numOfTrials = 192;
numOfUnits = 481;
numOfBins = 84;
meanOfFreq1 = zeros(numOfUnits, numOfBins);
meanOfFreq2 = zeros(numOfUnits, numOfBins);
meanOfFreq3 = zeros(numOfUnits, numOfBins);
meanOfFreq4 = zeros(numOfUnits, numOfBins);
meanOfFreq5 = zeros(numOfUnits, numOfBins);
meanOfFreq6 = zeros(numOfUnits, numOfBins);
s = zeros(numOfUnits, numOfTrials, numOfBins);
for k = 1 : numOfUnits
    r = randperm(numOfTrials);
    m1 = length(M(k).Cnd(1).TrialIdx);
    m2 = length(M(k).Cnd(1).TrialIdx) + length(M(k).Cnd(2).TrialIdx);
    m3 = length(M(k).Cnd(1).TrialIdx) + length(M(k).Cnd(2).TrialIdx) + length(M(k).Cnd(3).TrialIdx);
    m4 = length(M(k).Cnd(1).TrialIdx) + length(M(k).Cnd(2).TrialIdx) + length(M(k).Cnd(3).TrialIdx) + length(M(k).Cnd(4).TrialIdx);
    m5 = length(M(k).Cnd(1).TrialIdx) + length(M(k).Cnd(2).TrialIdx) + length(M(k).Cnd(3).TrialIdx) + length(M(k).Cnd(4).TrialIdx) + length(M(k).Cnd(5).TrialIdx);
    for i = 1 : numOfTrials
        temp = M(k).Trls{i, 1};
        for j = 1 : numOfBins
            s(k,i,j) =  length(temp(temp>= (j-1)*3.2/numOfBins-1.2 & temp<=(j)*3.2/numOfBins-1.2));
        end
        tmp(1:numOfBins) = s(k,i,:);
        if(ismember(i,r(1:m1)))
            meanOfFreq1(k, :) = meanOfFreq1(k,:) + tmp./length(M(k).Cnd(1).TrialIdx);
        elseif(ismember(i,r(m1+1:m2)))
            meanOfFreq2(k, :) = meanOfFreq2(k,:) + tmp./length(M(k).Cnd(2).TrialIdx);
        elseif(ismember(i,r(m2+1:m3)))
            meanOfFreq3(k, :) = meanOfFreq3(k,:) + tmp./length(M(k).Cnd(3).TrialIdx);
        elseif(ismember(i,r(m3+1:m4)))
            meanOfFreq4(k, :) = meanOfFreq4(k,:) + tmp./length(M(k).Cnd(4).TrialIdx);
        elseif(ismember(i,r(m4+1:m5)))
            meanOfFreq5(k, :) = meanOfFreq5(k,:) + tmp./length(M(k).Cnd(5).TrialIdx);
        elseif(ismember(i,r(m5+1:end)))
            meanOfFreq6(k, :) = meanOfFreq6(k,:) + tmp./length(M(k).Cnd(6).TrialIdx);
        end
    end
end
%% PCA
matrix = [meanOfFreq1(:, 33:72) meanOfFreq2(:, 33:72) meanOfFreq3(:, 33:72) meanOfFreq4(:, 33:72) meanOfFreq5(:, 33:72) meanOfFreq6(:, 33:72)];
%% Covariance matrix
covmat = zeros(6*40, 6*40);
for i = 1 :481
    covmat = covmat + ((matrix(i, :)- mean(matrix))')*(matrix(i, :)-mean(matrix));
end
covmat = covmat/480;
%% eigenvalues and eigenvectors
D = eig(covmat);
[V,~] = eig(covmat);
%% sort eigenvalues
[B,S] = sort(D, 'descend');
%%
close all
 figure
 for i = 1 : 6
temp3 = B(2)*V((i-1)*40+1:i*40, S(2));
temp4 = B(3)*V((i-1)*40+1:i*40, S(3));
plot(temp3,temp4)
hold on
 end
 xlabel("2nd principal component")
ylabel("3rd principal component")
title(["projection of shuffled Data to PC space", "84 bins for 3.2 sec"])
legend("CND1","CND2","CND3","CND4","CND5","CND6")






