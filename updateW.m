function [ W_new ] = updateW( W_old, edgeList )
% update W through adding edges
    W_new = W_old;
    for i = 1:length(edgeList)
        W_new(edgeList(i),edgeList(i)) = W_new(edgeList(i),edgeList(i)) + 1;
    end
end