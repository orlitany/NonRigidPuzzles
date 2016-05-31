function dy = diff_eta(x)
%     dy = 1-(tanh(2*x-1)).^2;
    dy = 3*(1- (tanh(6*x-3)).^2);
end