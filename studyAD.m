clear all; close all; clc;
cd('/Users/jyfeather/Google_Drive/Research/Manuscript/2016_BNLearn/codes');

%% load PET data 
load('/Volumes/mini/dataset/ADNI/PET.mat');
aal = csvread('/Volumes/mini/dataset/ADNI/AAL_YAN2.csv');

PET42 = PET(aal(:,1),:);

cohort = {'AD','MCI','NL'};
cohortNo = [49, 165, 232];

%% event variable via T2 chart
eventMat = zeros(232:42+3);                   % 42 metablism reduction events + 3 clinical events
eventMat(1:cohortNo(1),43) = 1;                % event 43 indicates AD
eventMat((cohortNo(1)+1):cohortNo(2),44) = 1;  % event 44 indicates MCI
eventMat((cohortNo(2)+1):cohortNo(3),45) = 1;  % event 45 indicates NL

for n = 1:42
    eventMat(:,n) = controlchart(PET42(n,:), 1);
end

%% learn orderings from event set
nEvals = 200;                        % a limit on the number of regression problems
discrete = 0;                        % Set to 1 for binary data
interv = 0;                          % Set to 0 for observational data
penalty = log(size(eventMat,1))/2;   % penalty according to BIC for learn DAG
potentialParents = ones(size(eventMat,2));
clamped = zeros(size(eventMat,1),size(eventMat,2));
numIter = 100;
nSamples = size(eventMat,1);

orderSet = [];
while numIter > 1
    numIter = numIter - 1;
    % get order
    [~,~,~,order,~,~] = OrderSearch(eventMat(randi(nSamples,1,nSamples),:),nEvals,0,penalty,discrete,clamped,potentialParents);
    orderSet = [orderSet; order];
end

%% compute the mean and covariance of ordering lambda0
mu0 = mean(orderSet, 1);
cov0 = cov(orderSet);
var0 = eye(size(cov0, 1));
lambda0 = pinv(cov0 / var0);

%% learn orderings using expert knowledge

%% draw ordering plots

%% draw 
