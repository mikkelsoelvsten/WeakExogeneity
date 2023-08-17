function [y,X] = genData(T,K,rho,a,S,model)
    % Error terms
    eps = randn(T+S,1);
    u = eps(S:end-1);
    for s=2:S
        u = u + eps(S+1-s:end-s)/s;
    end
    % Strictly exogenous part of regressors (orthogonalized)
    X = genRegressors(T,K,rho,model);
    % Weakly exogenous regressors    
    X(:,1) = X(:,1) + a * u;
    % Regression coefficients
    beta = zeros(K,1);
    % Outcomes
    y = X * beta + eps(S+1:end);
end

function [X] = genRegressors(T,K,rho,model)
    if isequal(model,'AR')
        % AR model for regressors
        mod = regARIMA('Intercept',0,'AR',{rho},'Var',1);
        % Generate regressors
        W = simulate(mod,T,'NumPaths',K,'U0',randn(T)/sqrt(1-rho^2));
    else
        % MA model for regressors
        mod = regARIMA('Intercept',0,'MA',{rho},'Var',1);
        % Generate regressors 
        W = simulate(mod,T,'NumPaths',K,'U0',randn(T)*sqrt(1+rho^2));
    end
    % Rotate regressors
    X = W / chol( W'*W/T );
end


