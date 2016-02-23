clear all; close all; clc;
cd('/Users/jyfeather/Google_Drive/Research/Manuscript/2016_BNLearn/codes');

%% setting need change before running experiments
numRun = 20;         % run 20 times to get mean and CI
data = 'asia';       % network from bayesian network repository
numExpert = 4;       % number of experts
sigma2 = 1;          % sigma square 

%% variables to record results
corr_SDP = [];      % numRun X numIter
var_SDP  = [];
corr_RD = [];
var_RD = [];

while numRun > 0
    clc;
    disp(numRun);
    numRun = numRun - 1;
    
    %% some parameters for observational data generation
    dagFunc = str2func(strcat('getDAG',data));
    [~, nodeNames] = dagFunc();
    numNodes = size(nodeNames, 2);
    nSamples = 2*numNodes;   % # of samples is twice # of nodes
     
    nEvals = 2500;           % a limit on the number of regression problems
    discrete = 0;            % Set to 1 for binary data
    interv = 0;              % Set to 0 for observational data

    [X,clamped,DAG,nodeNames] = sampleNetwork(data,nSamples,discrete,interv);

    %% some parameters for ordering finding
    penalty = log(nSamples)/2; % penalty according to BIC for learn DAG
    potentialParents = ones(size(X,2));

    %% initial order set
    [~,~,~,order,~,~] = OrderSearch(X,nEvals,0,penalty,discrete,clamped,potentialParents); orderSet = order;

    %% calculate order set
    numIter = 100;
    %[~,~,~,orderSet,~,~] = bootstrp(numIter, OrderSearch(X(randi(nSamples,1,nSamples),:),nEvals,0,penalty,discrete,clamped,potentialParents));
    while numIter > 1
        numIter = numIter - 1;
        % get order
        [~,~,~,order,~,~] = OrderSearch(X(randi(nSamples,1,nSamples),:),nEvals,0,penalty,discrete,clamped,potentialParents);
        orderSet = [orderSet; order];
    end

    %% compute the mean and covariance of ordering lambda0
    %orderSet = normr(orderSet);
    mu0 = mean(orderSet, 1);
    cov0 = cov(orderSet);
    var0 = eye(size(cov0, 1));
    lambda0 = pinv(cov0 / var0);
    % orderSet matrix: num of observations X num of nodes

    %% ground truth order phi, and arc incidence matrix B
    G = digraph(DAG, nodeNames); % construct DAG based on incidence matrix and nodes names
    phi = toposort(G);           % get ground truth order given DAG
    % DAG: row index = head, col index = tail
    % B: N rows(edges), n cols(nodes)
    n = size(DAG, 1);
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

    %% initial expert W0
    W0 = zeros(N);
    compare_pos = randi([1,N], 1, numExpert); % which edge is to be compared
    for i = 1:length(compare_pos)
        W0(compare_pos(i),compare_pos(i)) = 1;
    end

    %% generate expert scores y0
    y0 = mvnrnd(B*phi', sigma2*W0);

    %% update order by SDP
    numIter = 10;
    mu = [mu0; updateMu(mu0, W0, lambda0, y0, B)'];
    W = W0; lambda = lambda0;
    varPost = mean(diag(pinv(lambda)));
    %edgeCompared = [];
    for i=1:numIter
        % solve SDP
        edge2compare = SDP(W, lambda, B, n, N, numExpert);
        %edge2compare = randi([1,N],1,numExpert); % random selection
        %edgeCompared = [edgeCompared; edge2compare];

        % update W
        W = updateW(W, edge2compare);

        % generate y
        Winv = zeros(N);
        for j = 1:N
            if W(j,j)~=0
                Winv(j,j) = 1/W(j,j);
            end 
        end

        y = mvnrnd(B*phi', sigma2*Winv);

        % update mu, which is the updated order
        mu = [mu; updateMu(mu(end,:), W, lambda, y, B)'];

        % update lambda and posterior variance
        lambda = B'*W*B + lambda;
        varPost = [varPost mean(diag(pinv(lambda)))];
    end
    
    corrSet = [];
    for i = 1:numIter
        corrSet = [corrSet corr(phi', mu(i,:)')];
    end
    
    corr_SDP = [corr_SDP; corrSet];
    var_SDP = [var_SDP; varPost(1:numIter)];

    %% update order by random selection
    numIter = 10;
    mu = [mu0; updateMu(mu0, W0, lambda0, y0, B)'];
    W = W0; lambda = lambda0;
    varPost = mean(diag(pinv(lambda)));
    %edgeCompared = [];
    for i=1:numIter
        edge2compare = randi([1,N],1,numExpert); % random selection
        %edgeCompared = [edgeCompared; edge2compare];

        % update W
        W = updateW(W, edge2compare);

        % generate y
        Winv = zeros(N);
        for j = 1:N
            if W(j,j)~=0
                Winv(j,j) = 1/W(j,j);
            end 
        end

        y = mvnrnd(B*phi', sigma2*Winv);

        % update mu, which is the updated order
        mu = [mu; updateMu(mu(end,:), W, lambda, y, B)'];

        % update lambda and posterior variance
        lambda = B'*W*B + lambda;
        varPost = [varPost mean(diag(pinv(lambda)))];
    end
    
    corrSet = [];
    for i = 1:numIter
        corrSet = [corrSet corr(phi', mu(i,:)')];
    end
    
    corr_RD = [corr_RD; corrSet];
    var_RD = [var_RD; varPost(1:numIter)];
end

%% save result to csv for visualization
csvwrite(['./result/', data, '_', int2str(numExpert), '_', int2str(sigma2), '_corr_SDP.csv'], corr_SDP);
csvwrite(['./result/', data, '_', int2str(numExpert), '_', int2str(sigma2), '_var_SDP.csv'], var_SDP);
csvwrite(['./result/', data, '_', int2str(numExpert), '_', int2str(sigma2), '_corr_RD.csv'], corr_RD);
csvwrite(['./result/', data, '_', int2str(numExpert), '_', int2str(sigma2), '_var_RD.csv'], var_RD);
