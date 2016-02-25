clear all; close all; clc;
cd('/Users/jyfeather/Google_Drive/Research/Manuscript/2016_BNLearn/codes');

%% load KPI data 
KPI = csvread('../data/KPI.csv');

KPInames = {'pm_ermin','pm_hsmbn','pm_hsvcn','pm_hsvgn','sp_hsdmn','sp_hsesn','sp_hsifn','sp_hsmfn',...
'sp_hsmgn','sp_hspcn','sp_hsrsn02','sp_hsrsn03','sp_hsrsn04','sp_hsrsn05',...
'sp_hsrsn06','sp_hsrsn07','sp_hsrsn08','sp_hsrsn09','sp_hsrsn10',...
'sp_hstmn','sp_hstwn01','sp_psltn','sp_pssdn','sp_psspn','sp_pstrn',...
'hr_hrcrn','hr_hrtnn','hr_hrtrn02','hr_hrtrn13','hr_hrtrn04','hr_irqin'};

%% learn orderings from event set
nEvals = 200;                        % a limit on the number of regression problems
discrete = 0;                        % Set to 1 for binary data
interv = 0;                          % Set to 0 for observational data
penalty = log(size(KPI,1))/2;   % penalty according to BIC for learn DAG
potentialParents = ones(size(KPI,2));
clamped = zeros(size(KPI,1),size(KPI,2));
numIter = 100;
nSamples = size(KPI,1);

orderSet = [];
while numIter > 1
    numIter = numIter - 1;
    % get order
    [adj,~,~,order,~,~] = OrderSearch(KPI(randi(nSamples,1,nSamples),:),nEvals,0,penalty,discrete,clamped,potentialParents);
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
numExpert = 50;
edge2compare = SDP(W0, lambda0, B, n, N, numExpert);

%% write edge2compare into csv
compareEdge = {};
for iter = 1:numExpert
    noHead = find(B(edge2compare(iter),:)==1);
    noTail = find(B(edge2compare(iter),:)==-1);
    compareEdge = [compareEdge; {KPInames(noHead), KPInames(noTail)}];
end

%% draw ordering plots
[~, ordering0] = sort(mu0);
[~, orderingFinal] = sort(ordering0);

%% draw 
