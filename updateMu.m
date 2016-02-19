function [ mu_new ] = updateMu( mu_old, W, lambda, y, B )
% update mu
    mu_new = pinv(B'*W*B+lambda)*(B'*W*y'+lambda*mu_old');
end




