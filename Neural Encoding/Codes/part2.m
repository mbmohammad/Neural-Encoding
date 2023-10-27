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
numOfBins = 84;
meanOfFreq1 = zeros(numOfUnits, numOfBins);
meanOfFreq2 = zeros(numOfUnits, numOfBins);
meanOfFreq3 = zeros(numOfUnits, numOfBins);
meanOfFreq4 = zeros(numOfUnits, numOfBins);
meanOfFreq5 = zeros(numOfUnits, numOfBins);
meanOfFreq6 = zeros(numOfUnits, numOfBins);
s = zeros(numOfUnits, numOfTrials, numOfBins);
for k = 1 : numOfUnits
    for i = 1 : numOfTrials
        temp = M(k).Trls{i, 1};
        for j = 1 : numOfBins
            s(k,i,j) =  length(temp(temp>= (j-1)*3.2/numOfBins-1.2 & temp<=(j)*3.2/numOfBins-1.2));
        end
        tmp(1:numOfBins) = s(k,i,:);
        if(ismember(i,M(k).Cnd(1).TrialIdx))
            meanOfFreq1(k, :) = meanOfFreq1(k,:) + tmp./length(M(k).Cnd(1).TrialIdx);
        elseif(ismember(i,M(k).Cnd(2).TrialIdx))
            meanOfFreq2(k, :) = meanOfFreq2(k,:) + tmp./length(M(k).Cnd(2).TrialIdx);
        elseif(ismember(i,M(k).Cnd(3).TrialIdx))
            meanOfFreq3(k, :) = meanOfFreq3(k,:) + tmp./length(M(k).Cnd(3).TrialIdx);
        elseif(ismember(i,M(k).Cnd(4).TrialIdx))
            meanOfFreq4(k, :) = meanOfFreq4(k,:) + tmp./length(M(k).Cnd(4).TrialIdx);
        elseif(ismember(i,M(k).Cnd(5).TrialIdx))
            meanOfFreq5(k, :) = meanOfFreq5(k,:) + tmp./length(M(k).Cnd(5).TrialIdx);
        elseif(ismember(i,M(k).Cnd(6).TrialIdx))
            meanOfFreq6(k, :) = meanOfFreq6(k,:) + tmp./length(M(k).Cnd(6).TrialIdx);
        end
    end
end
%% GLM analysis first part(set X)
X1 = zeros(numOfUnits,numOfTrials);
X2 = zeros(numOfUnits,numOfTrials);
X3 = zeros(numOfUnits,numOfTrials);
X4 = zeros(numOfUnits,numOfTrials);
X5 = zeros(numOfUnits,numOfTrials);
X6 = zeros(numOfUnits,numOfTrials);
for k = 1 : numOfUnits
    X1(k, M(k).Cnd(1).TrialIdx) = 1;
    X2(k, M(k).Cnd(2).TrialIdx) = 1;
    X3(k, M(k).Cnd(3).TrialIdx) = 1;
    X4(k, M(k).Cnd(4).TrialIdx) = 1;
    X5(k, M(k).Cnd(5).TrialIdx) = 1;
    X6(k, M(k).Cnd(6).TrialIdx) = 1;
end
%% GLM analysis 2nd part(set Y)
Y = zeros(numOfUnits, numOfTrials, 15);
for k = 1 : numOfUnits
    for i = 1 : numOfTrials
        Y(k, i, :) = s(k,i,12:26);
    end
end
%% GLM analysis 3rd part(regress)
coeffs = zeros(6, numOfUnits, 15);
Stats = zeros(4, numOfUnits, 15);
for k = 1 : numOfUnits
    for i = 1 : 15
        X = [X1(k,:)' X2(k,:)' X3(k,:)' X4(k,:)' X5(k,:)' X6(k,:)'];
        Yreg = Y(k, :, i)';
        [b,bint,r,rint,stats] = regress(Yreg,X);
        coeffs(:,k, i) = b;
        Stats(:, k, i) = stats;
    end
end
%% sum of P-Values
sumStats = sum(Stats, 3);
%% sort stats
sumStatsT = sumStats';
[sortedStats, idxx] = sort(sumStatsT(:, 3));
%%



%% TODO
close all
xaxis = -0.1: 0.1 :1.3;
for k = 457 : 457
    iiii(:,:) = coeffs(:, k,:);
    figure
    plot(xaxis, iiii')
    xlabel("time compared to cue onset(sec)")
    ylabel("coefficient")
    title(sprintf("coefficients of GLM analysis for neuron number %d", k))
    legend("CND1","CND2","CND3","CND4","CND5","CND6")
end

%% PCA
% matrix = [meanOfFreq1(:, 33:72) meanOfFreq2(:, 33:72) meanOfFreq3(:, 33:72) meanOfFreq4(:, 33:72) meanOfFreq5(:, 33:72) meanOfFreq6(:, 33:72)];
matrix = meanOfFreq6(:, 33:72);
%% Covariance matrix
covmat = zeros(40, 40);
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
hold on
temp3 = B(3)*V(:, S(3));
temp4 = B(2)*V(:, S(2));
temp5 = B(4)*V(:, S(4));
plot3(temp3,temp4, temp5)
%%

legend("CND1","CND2","CND3","CND4","CND5","CND6")
zlabel("4th principal component")
ylabel("2nd principal component")
xlabel("3rd principal component")
title(["projection of Data to PC space", "84 bins for 3.2 sec"])



%%
close all
 figure
 for i = 1 : 6
temp3 = B(10)*V((i-1)*40+1:i*40, S(10));
temp4 = B(6)*V((i-1)*40+1:i*40, S(6));
plot(temp3,temp4)
hold on
 end
 xlabel("4th principal component")
ylabel("3rd principal component")
title(["projection of Data to PC space", "84 bins for 3.2 sec"])
legend("CND1","CND2","CND3","CND4","CND5","CND6")

%% 3d plot
figure
for i = 1 :6
temp3 = V((i-1)*40+1:i*40, S(1));
temp4 = V((i-1)*40+1:i*40, S(2));
temp5 = V((i-1)*40+1:i*40, S(3));
plot3(temp3, temp4, temp5)
hold on
end
xlabel("1st principal component")
ylabel("2nd principal component")
zlabel("3rd principal component")
title(["projection of Data to PC space", "84 bins for 3.2 sec"])
legend("CND1","CND2","CND3","CND4","CND5","CND6")






%%
close all
 figure
temp3 = B(3)*V(1:40, S(3));
temp4 = B(2)*V(1:40, S(2));
plot(temp3,temp4)
hold on
temp3 = B(3)*V(41:80, S(3));
temp4 = B(2)*V(41:80, S(2));
plot(temp3,temp4)
hold on
temp3 = B(3)*V(81:120, S(3));
temp4 = B(2)*V(81:120, S(2));
plot(temp3,temp4)
hold on
temp3 = B(3)*V(121:160, S(3));
temp4 = B(2)*V(121:160, S(2));
plot(temp3,temp4)
hold on
temp3 = B(3)*V(161:200, S(3));
temp4 = B(2)*V(161:200, S(2));
plot(temp3,temp4)
hold on
temp3 = B(3)*V(201:240, S(3));
temp4 = B(2)*V(201:240, S(2));
plot(temp3,temp4)
xlabel("3rd principal component")
ylabel("2nd principal component")
title(["projection of Data to PC space", "84 bins for 3.2 sec"])
legend("CND1","CND2","CND3","CND4","CND5","CND6")
%%
