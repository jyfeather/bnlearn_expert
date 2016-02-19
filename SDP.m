function [ edge2compare ] = SDP( W, lambda, B, n, N, numExpert )
% solve SDP to get edge list to be compared in the next time
% input:
%   - W:                matrix to record compare times
%   - lambda:           inverse covariance matrix
%   - B:                arc incidence matrix
%   - n:                num of nodes
%   - N:                num of potential edges
%
% output:
%   - edge2compare:     the edge list

epsilon = numExpert;         % number of edges need to compare
lambda = (lambda'+lambda)/2; % lambda might be not symmetric because of pinv

% solve SDP
cvx_begin sdp
    % variable declare
    variable s;
    variable x(N) nonnegative;
    expression z;
    z = zeros(n);
    for i = 1:N
        z = z + x(i)*(B(i,:)'*B(i,:));
    end
    
    % formulation
    maximize(s);
    subject to
        s*(eye(n)-ones(n)/n) <= B'*W*B + lambda + z;
        %s*(eye(n)-ones(n)/n) <= L + z;
        sum(x) <= epsilon;
        %x(edgeCompared) == 1;
        x <= 1;
cvx_end

[~, edge2compare] = sort(x, 'descend');
edge2compare = edge2compare(1:epsilon);

end % end of function