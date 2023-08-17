function [g,Bias] = tr(D,M,DP,T,S,gam1,type)
% Function to compute traces entering the bias

    if strcmp(type,'future')
        gam2 = gam1;
        gam1 = zeros(length(gam1),1);
    elseif strcmp(type,'both')
        gam2 = gam1;
    elseif strcmp(type,'free')
        gam2 = gam1(1+length(gam1)/2:end);
        gam1 = gam1(1:length(gam1)/2);
    else
        gam2 = zeros(length(gam1),1);
    end

    A = eye(T) - gam1(1) * D{1};
    B = eye(T) + gam2(1) * D{1} - ( gam1(1) + gam2(1) ) * DP{1};
    for s = 2:length(gam1)
        A = A - gam1(s) * D{s};
        B = B + gam2(s) * D{s} - ( gam1(s) + gam2(s) ) * DP{s};
    end
    Bias = M * (B \ A);
    g = sum( diag( Bias, -1 ) );
    for s = 2:S
        g = [ g ; sum( diag( Bias, -s ) ) ];
    end
end
