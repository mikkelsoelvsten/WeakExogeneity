function [OLS,IV,trace,BbOLS,BbIV] = reg_func(y,X,S,type,coef,Sb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explanations:
% This function computes and returns (i) OLS estimates and (ii) IV estimates, 
% which are constructed to alleviate bias from many weakly exogeneous regressors.
% Additional outputs are (i) 'tr', which computes the traces that determine
% the potential for bias in OLS.
%
% The inputs 'y' and 'X' are the vector of outcomes and matrix of regressors, respectively.
% The input 'S' specifies the number of periods for which the IV estimator
% should account for a failure in strict exogeneity
% The input 'type' can be 'lead' or 'future'. When 'type' is 'lead' the IV
% estimator uses the regressor and its S leads to construct a technical
% instrument. When 'type' is 'future' the IV estimator uses the regressors
% and all its leads to construct an instrument.
% The input 'coef' specifies which of the coefficients to return.
% 'type' and 'coef' are optional with 'lead' and all coefficients as defaults.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BbOLS and BBIV computes the bias bounds for OLS and IV - these bounds are computed as a fraction of
% standard deviation and for up to 'Sb' periods of feedback.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set optional parameters if not provided
if nargin <= 3
    type = 'lead';
end
if nargin <= 4
    coef = 1:size(X,2);
end
if nargin <= 5
    Sb = 20;
end


%% Preliminaries
T = length(X);
R = (X'*X)\X';
P = X * R;
M = eye(T) - P;
for s=1:Sb
D{s} = spdiags(ones(T,1),-s,T,T);
end
for s=1:S
DP{s}= D{s} * P;
Dy{s}= D{s} * y;
end


%% The observable parts of the bias in OLS
for s=1:Sb
    trace{s} = sum(diag(M,-s));
end


%% OLS
OLS = R*y;
OLS = OLS(coef);


%% IV
if S>1 || strcmp(type,'free')
    % Impose the 2^S inequality constraints that the 1-norm of gamma <= 1
    options = optimoptions('fmincon','Display','off',...
        'MaxFunctionEvaluations',3e4,'OutputFcn',@outfun);
    A = [1,-1];
    for s = 1:S-1
        A = combvec(A,[1,-1]);
    end
    b = .99*ones(2^S,1);

    % Compute gamma
    fun1 = @(x) sum(tr(D,M,DP,T,S,x,type).^2);
    init = zeros(S,1);

    if strcmp(type,'free')
        A = blkdiag(A,A);
        b = [b;b];
        init = [init;init];
    end
    
    
    g = fmincon( fun1, init, A',b,[],[],[],[],[],options);

    if strcmp(type,'future')
        g1 = zeros(length(g),1);
        g2 = g;
    elseif strcmp(type,'both')
        g1 = g;
        g2 = g;
    elseif strcmp(type,'free')
        g1 = g(1:length(g)/2);
        g2 = g(1+length(g)/2:end);
    else
        g1 = g;
        g2 = zeros(length(g),1);
    end

    A = eye(T);
    B = y;
    for s = 1:S
        A = A + g2(s) * D{s} - ( g1(s) + g2(s) ) * DP{s};
        B = B - g1(s) * Dy{s};
    end

    % Transformed outcomes
    y1 = A \ B;
else
    % Compute gamma
    fun1 = @(x) tr(D,M,DP,T,S,x,type)^2;
    g = fminbnd( fun1, -.99,.99);

    if strcmp(type,'future')
        g1 = zeros(length(g),1);
        g2 = g;
    elseif strcmp(type,'both')
        g1 = g;
        g2 = g;
    else
        g1 = g;
        g2 = zeros(length(g),1);
    end

    % Transformed outcomes
    y1 = ( eye(T) + g2 * D{1} - ( g1 + g2 ) * DP{1} ) \ ( y - g1 * Dy{1} );

end

IV = R*y1;
IV = IV(coef);


%% Bias bounds
if nargout >= 4

%% Preliminaries
for s=1:Sb
    DpP{s}= D{s}' * P;
end

%% Bias bound for OLS

% The observable parts of the bias in OLS
traces = cell2mat(trace);

% Computing (the lower bound on) the variance matrix
A = {};
for s=1:Sb
    A{s} = (D{s} + D{s}' - DpP{s} - DpP{s}')/2;
end
V = zeros(Sb,Sb);
for s = 1:Sb
    for t = s:Sb
        V(s,t) = 2*sum( A{s} .* A{t}, 'all');
        V(t,s) = V(s,t);
    end
end

Bb = traces / chol( V );
BbOLS = cumsum( Bb.^2 ).^(.5)';

%% Bias bound for IV

% The observable parts of the bias in IV
[traces,Mg] = tr(D,M,DP,T,Sb,g,type);

% Computing (the lower bound on) the variance matrix
for s=1:Sb
    DMg{s} = D{s}' * Mg;
end
for s=1:Sb
    A{s} = ( DMg{s} + DMg{s}' )/2;
end
V = zeros(Sb,Sb);
for s = 1:Sb
    for t = s:Sb
        V(s,t) = 2*sum( A{s} .* A{t}, 'all');
        V(t,s) = V(s,t);
    end
end
Bb = traces' / chol( V );
BbIV = cumsum( Bb.^2 ).^(.5)';

end


end

%% Functions

% To control the optimization procedure
function stop = outfun(~,optimValues,~)
stop = false;
% Check if objective function is less than 1e-4.
if optimValues.fval < 1e-4^2
    stop = true;
end
end
