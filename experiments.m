clear all; close all; clc;
cd('/Users/jyfeather/Google_Drive/Research/Manuscript/2016_BNLearn/codes');

%% all possible combinations
dataNameSet = {'earthquake', 'asia', 'child', 'insurance', 'mildew', 'alarm', 'barley', 'hailfinder'};
dataNumSet = 1:size(dataNameSet,2);
numExpertSet = [1 2 4];       
sigma2Set = [1 2 4];       
combinations = combvec(dataNumSet, numExpertSet, sigma2Set); % all possible combinations

numComb = 30;
while numComb <= size(combinations,2)
    try
    clc;
    disp(numComb);
        
    %% prepare some settings
    eachComb = combinations(:,numComb);
    data = dataNameSet{eachComb(1)};
    numExpert = eachComb(2);
    sigma2 = eachComb(3);
    numRun = 20;         
    
    %% variables to record results
    corr_SDP = [];      % numRun X numIter
    var_SDP  = [];
    corr_RD = [];
    var_RD = [];

    while numRun > 0
        numRun = numRun - 1;
        disp(numRun);
        
        %% some parameters for observational data generation
        dagFunc = str2func(strcat('getDAG',data));
        [~, nodeNames] = dagFunc();
        numNodes = size(nodeNames, 2);
        nSamples = 2*numNodes;   % # of samples is twice # of nodes

        nEvals = 20;           % a limit on the number of regression problems
        discrete = 0;            % Set to 1 for binary data
        interv = 0;              % Set to 0 for observational data

        [X,clamped,DAG,nodeNames] = sampleNetwork(data,nSamples,discrete,interv);

        %% some parameters for ordering finding
        penalty = log(nSamples)/2; % penalty according to BIC for learn DAG
        potentialParents = ones(size(X,2));

        %% initial order set
        [~,~,~,order,~,~] = OrderSearch(X,nEvals,0,penalty,discrete,clamped,potentialParents); 
        orderSet = order;

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

        %% update order by SDP
        numIter = 10;
        W0 = zeros(N);
        mu = mu0;
        W = W0; 
        lambda = lambda0;
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
        W0 = zeros(N);
        mu = mu0;
        W = W0; 
        lambda = lambda0;
        varPost = mean(diag(pinv(lambda)));
        for i=1:numIter
            edge2compare = randi([1,N],1,numExpert); % random selection

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
    
    numComb = numComb + 1;
    catch
    % do nothing
    end
end
