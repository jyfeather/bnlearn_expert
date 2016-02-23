function [ res ] = controlchart( vec, numSigma )
    x_bar = mean(vec);
    x_sd = std(vec);
    x_min = x_bar-numSigma*x_sd;
    x_max = x_bar+numSigma*x_sd;
    
    res_min = (vec < x_min);
    res_max = (vec > x_max);
    res = res_min | res_max;
end

