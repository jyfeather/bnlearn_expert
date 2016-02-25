clear all; close all; clc;
cd('/Users/jyfeather/Google_Drive/Research/Manuscript/2016_BNLearn/codes');

%% load PET data 
load('/Volumes/mini/dataset/ADNI/PET.mat');
aal = csvread('/Volumes/mini/dataset/ADNI/AAL_YAN2.csv');
ROInames = {'Frontal_Sup','Frontal_Mid',...
'Frontal_Sup_Medial','Frontal_Mid_Orb','Rectus','Cingulum_Ant','Parietal_Sup',...
'Parietal_Inf','Precuneus','Cingulum_Post','Occipital_Sup','Occipital_Mid',...
'Occipital_Inf','Temporal_Sup','Temporal_Pole_Sup','Temporal_Mid',...
'Temporal_Pole_Mid','Temporal_Inf','Fusiform','Hippocampus','ParaHippocampal'};

PET42 = PET(aal(:,1),:);

cohort = {'AD','MCI','NC'};
cohortNo = [49, 165, 232];

%% event variable via T2 chart
eventMat0 = zeros(232,42);                     % 42 metablism reduction events
for n = 1:42
    eventMat0(:,n) = controlchart(PET42(n,:), 1);
end

% measurements on same regions on left and right side is converted into one
eventMat = zeros(232,21);
for n = 1:21
    eventMat(:,n) = eventMat0(:,2*n-1) | eventMat0(:,2*n);
end

%  21 + 3 clinical events
eventMat(1:cohortNo(1),22) = 1;                % event 22 indicates AD
eventMat((cohortNo(1)+1):cohortNo(2),23) = 1;  % event 23 indicates MCI
eventMat((cohortNo(2)+1):cohortNo(3),24) = 1;  % event 24 indicates NL

eventNames = [ROInames, 'AD', 'MCI', 'NC'];

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
    [adj,~,~,order,~,~] = OrderSearch(eventMat(randi(nSamples,1,nSamples),:),nEvals,0,penalty,discrete,clamped,potentialParents);
    orderSet = [orderSet; order];
end

%% compute the mean and covariance of ordering lambda0
mu0 = mean(orderSet, 1);
cov0 = cov(orderSet);
var0 = eye(size(cov0, 1));
lambda0 = pinv(cov0 / var0);

%% learn orderings using expert knowledge
%% ground truth order phi, and arc incidence matrix B
% B: N rows(edges), n cols(nodes)
n = size(mu0,2);
N = nchoosek(n, 2);
B = zeros(N, n);
B_pos = 1;
for i = 1:(n-1)
    for j = (i+1):n
        B(B_pos,i) = 1;  % i = head
        B(B_pos,j) = -1; % j = tail
        B_pos = B_pos + 1;
    end 
end
W0 = zeros(N);
numExpert = 20;
edge2compare = SDP(W0, lambda0, B, n, N, numExpert);

%% write edge2compare into csv
compareEdge = {};
for iter = 1:numExpert
    noHead = find(B(edge2compare(iter),:)==1);
    noTail = find(B(edge2compare(iter),:)==-1);
    compareEdge = [compareEdge; {eventNames(noHead), eventNames(noTail)}];
end

%% draw ordering plots
[~, ordering0] = sort(mu0);
[~, orderingFinal] = sort(ordering0);

%% draw 
