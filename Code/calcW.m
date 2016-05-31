function [ W ] = calcW( k,est_rank )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
W = zeros(k);
for i=1:k
    for j=1:k
        slope = est_rank/k;
        direction = [1 slope];
        direction = direction./norm(direction);
        W(i,j) = exp(-0.03*sqrt(i.^2 + j.^2))*norm(cross([direction 0], [i,j, 0]-[1 1 0]));        
    end
end

end

