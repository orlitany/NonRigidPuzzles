function [D, matches] = run_icp(M, N, r, C_init, max_iters)

fprintf('Running ICP...\n');

C_init = C_init(:,1:r);

X = N.evecs(:,1:r)' * N.S;
Y = M.evecs' * M.S;
tree = ann('init', Y);

matches = ann('search', tree, C_init*X, 1)';

err = sum( sqrt(sum((C_init*X - Y(:,matches)).^2)) );
err = err / (size(X,2)*size(C_init,1));

fprintf('(0) MSE: %.2e\n', err);

if max_iters == 0
    ann('deinit', tree);
    D = C_init;
    return
end

% Start iterations

D_prev = C_init;
err_prev = err;

for i=1:max_iters
    
    [U,~,V] = svd(X * Y(:,matches)');
    D = U * V(:,1:r)';
    D = D';
    
    matches = ann('search', tree, D*X, 1)';

    err = sum( sqrt(sum((D*X - Y(:,matches)).^2)) );
    err = err / (size(X,2)*size(C_init,1));
    
    fprintf('(%d) MSE: %.2e\n', i, err);
    
    if err > err_prev
        fprintf('Local optimum reached.\n');
        D = D_prev;
        break;
    end
    
    if (err_prev - err) < 5e-6
        fprintf('Local optimum reached.\n');
        break;
    end
    
    err_prev = err;
    D_prev = D;

end

ann('deinit', tree);
