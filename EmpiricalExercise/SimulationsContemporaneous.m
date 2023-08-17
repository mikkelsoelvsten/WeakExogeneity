clearvars

%% Load data set of 204 observations and 108 variables
% We will be using no lags in the analysis but will work with 200 total
% observations
[analysisData,analysisData_names] = LoadAnalysisData;
T = size(analysisData,1)-4;
Ktotal = size(analysisData,2);

%% Setting up simulations for randomly choosing one outcome and K1 regressors to use contemporaneously
rep = 1000;
dgps = 100;
Ks = [5 15 25 35 45 55 65 75 85 100];
Knum = length(Ks);
Ksi = Ks(1):0.1:Ks(Knum);


IRF = zeros(dgps,Knum);
IRFiv = zeros(dgps,Knum);
IRFse = zeros(dgps,Knum);
IRFivse = zeros(dgps,Knum);
IRFdiffse = zeros(dgps,Knum);
COEF = zeros(dgps,Ks(Knum),Knum);
COEFiv = zeros(dgps,Ks(Knum),Knum);
SE = zeros(dgps,Ks(Knum),Knum);
SEiv = zeros(dgps,Ks(Knum),Knum);
SEDIFF = zeros(dgps,Ks(Knum),Knum);
Ftest = zeros(dgps,Knum);
ESD = zeros(dgps,Knum);
TRACE = zeros(dgps,Knum);
BIASFORMULA = zeros(dgps,Knum);

MEANs = zeros(dgps,6,Knum);
STDs = zeros(dgps,6,Knum);
ESDs = zeros(dgps,2,Knum);
Ts = zeros(dgps,2,Knum);
T1s = zeros(dgps,2,Knum);

%% 

tic
for k=1:Knum
    K = Ks(k)
    parfor i=1:dgps
        % Fix seed for random number generation
        rng(k*100+i);

        % Randomly select variables for regression
        v = rand(Ktotal,1);
        sorted_v = sort(unique(v), 'descend');
        subsety = (v >= sorted_v(1));
        subsetX = (v < sorted_v(1)) & (v >= sorted_v( K ));
        yt = analysisData(5:end,subsety);
        Xt = [ones(T,1),analysisData(5:end,subsetX)];
        
        % Compute parameters for the DGP
        A = Xt'*Xt;
        b = A\Xt'*yt;
        [~,biv,~,~,~,~,~,~,~,~,~,VLS,VIV,VDIFF] ...
            = reg_func(yt,Xt,1,'lead');
        e = yt - Xt * b;
        s = sqrt(e'*e/(T-K));
        a = Xt(2:end,:)'*e(1:end-1)/(e'*e);
        a(1) = 0;
        Bparm = Xt(2:end,:)'*Xt(1:end-1,:) / (Xt'*Xt);

        COEF(i,:,k) = [b; zeros(Ks(Knum)-K,1)];
        COEFiv(i,:,k) = [biv; zeros(Ks(Knum)-K,1)];
        SE(i,:,k) = [sqrt(diag(VLS)); zeros(Ks(Knum)-K,1)];
        SEiv(i,:,k) = [sqrt(diag(VIV)); zeros(Ks(Knum)-K,1)];
        SEDIFF(i,:,k) = [sqrt(diag(VDIFF)); zeros(Ks(Knum)-K,1)];
        Ftest(i,k) = (b-biv)'*(pinv(VDIFF)*(b-biv))/rank(VDIFF);
        IRF(i,k) = a'*b;
        IRFse(i,k) = sqrt(a'*VLS*a);
        IRFiv(i,k) = a'*biv;
        IRFivse(i,k) = sqrt(a'*VIV*a);
        IRFdiffse(i,k) = sqrt(a'*VDIFF*a);
        ESD(i,k) = s;
        TRACE(i,k) = - sum(diag(Bparm));
        BIASFORMULA(i,k) = s^2 * a'*((A+s^2*(T-K)*a*a')\a) * TRACE(i,k);

        % Conduct `rep' simulations for this design
        estimates = zeros(rep,6);
        stds = zeros(rep,2);
        ses = zeros(rep,2);
        for j=1:rep
            % Generate data from the specified DGP
            eps = s*randn(T+5,1);
            X = [ones(T,1),Xt(:,2:end) + eps(5:end-1) * a(2:end)'];
            y = X * b + eps(6:end);
            % Estimate IRF and error variance estimates
            [ols,iv,varLS,varIV,irf_lp,irf_LS,irf_IV,irf_IV2,~,~,~,covarLS,covarIV] ...
                = reg_func(y,X,1,'lead');
            estimates(j,:) = [ols'*a, iv'*a, irf_lp, irf_LS, irf_IV, irf_IV2];
            stds(j,:) = [sqrt(varLS) sqrt(varIV)];
            ses(j,:) = [sqrt(a'*covarLS*a) sqrt(a'*covarIV*a)];
        end

        % Mean IRFs and their standard deviations
        MEANs(i,:,k) = mean(estimates);
        STDs(i,:,k) = std(estimates);
        % Mean of variance estimates
        ESDs(i,:,k) = mean(stds);
        Ts(i,:,k) = mean( abs( (estimates(:,1:2)-IRF(i,k))./ses ) > 1.96 );
        T1s(i,:,k) = mean(  (estimates(:,1:2)-IRF(i,k))./ses < -1.64 );
    end
end

toc

workspaceName = 'SimulationsContemporaneous.mat';
save(fullfile('SavedResults',workspaceName));

%%

sim.IRF = IRF;
sim.ESD = ESD;
sim.TRACE = TRACE;
sim.BIASFORMULA = BIASFORMULA;

sim.MEANs = MEANs;
sim.STDs = STDs;
sim.ESDs = ESDs;
for k=1:Knum
    sim.Bias(:,:,k) = sim.MEANs(:,:,k) - sim.IRF(:,k);
end

%% Find the DGP of the 10th and 90th percentile OLS bias
absbias = abs(sim.Bias(:,1,:));
A = quantile(absbias,.105);
B10 = absbias <= A & absbias >= A;

A = quantile(absbias,.905);
B90 = absbias <= A & absbias >= A;


%% 
f1=figure('Renderer', 'painters', 'Position', [10 10 700 300]);

absiv = abs(sim.Bias(:,2,:));
stdls = sim.STDs(:,1,:);
stdiv = sim.STDs(:,2,:);

A1 = polyval( polyfit(Ks,absbias(B10),5), Ksi);
A2 = polyval( polyfit(Ks,absiv(B10),5), Ksi);
B1 = polyval( polyfit(Ks,stdls(B10),5), Ksi);
B2 = polyval( polyfit(Ks,stdiv(B10),5), Ksi);
C = polyval( polyfit(Ks,abs( sim.TRACE(B10) )/T,5), Ksi);

s1=subplot(2,2,1);
p1=plot(Ksi,A1,'r',Ksi,B1,'--r',Ksi,A2,'b',Ksi,B2,'--b',Ksi,C,':k')
ylim([-.01 .15])
xticks([5 25 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)
title("10th percentile experiment")

A1 = polyval( polyfit(Ks,absbias(B90),5), Ksi);
A2 = polyval( polyfit(Ks,absiv(B90),5), Ksi);
B1 = polyval( polyfit(Ks,stdls(B90),5), Ksi);
B2 = polyval( polyfit(Ks,stdiv(B90),5), Ksi);
C = polyval( polyfit(Ks,abs( sim.TRACE(B90) )/T,5), Ksi);

s2=subplot(2,2,2);
p2=plot(Ksi,A1,'r',Ksi,B1,'--r',Ksi,A2,'b',Ksi,B2,'--b',Ksi,C,':k')
ylim([-.01 .15])
xticks([5 25 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)
title("90th percentile experiment")

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Bias (OLS)','Stddev (OLS)','Bias (IV)','Stddev (IV)','Trace/$T$','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
axis(hLegend,'off');
box(leg,'off');
set(leg,'Position',posLegend);

% Adjusting positions of first two plots
poss1 = get(s1,'Position');
poss1(2) = poss1(2) - .37;
poss1(4) = poss1(4) + .37;
set(s1,'Position',poss1);

poss2 = get(s2,'Position');
poss2(2) = poss2(2) - .37;
poss2(4) = poss2(4) + .37;
set(s2,'Position',poss2);

%%
f2=figure('Renderer', 'painters', 'Position', [10 10 700 300]);

Tsls = Ts(:,1,:);
Tsiv = Ts(:,2,:);

A1 = polyval( polyfit(Ks,Tsls(B10),4), Ksi);
A2 = polyval( polyfit(Ks,Tsiv(B10),4), Ksi);

%A1 = Tsls(B10);
%A2 = Tsiv(B10);

s1=subplot(2,2,1);
p1=plot(Ksi,A1,'r',Ksi,A2,'b')
ylim([.0 .35])
xticks([5 25 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)
title("10th percentile experiment")

A1 = polyval( polyfit(Ks,Tsls(B90),4), Ksi);
A2 = polyval( polyfit(Ks,Tsiv(B90),4), Ksi);

%A1 = Tsls(B90);
%A2 = Tsiv(B90);

s2=subplot(2,2,2);
p2=plot(Ksi,A1,'r',Ksi,A2,'b')
ylim([.0 .35])
xticks([5 25 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)
title("90th percentile experiment")

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Size (OLS)','Size (IV)','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
axis(hLegend,'off');
box(leg,'off');
set(leg,'Position',posLegend);

% Adjusting positions of first two plots
poss1 = get(s1,'Position');
poss1(2) = poss1(2) - .37;
poss1(4) = poss1(4) + .37;
set(s1,'Position',poss1);

poss2 = get(s2,'Position');
poss2(2) = poss2(2) - .37;
poss2(4) = poss2(4) + .37;
set(s2,'Position',poss2);


%% Find the DGP of the 10th and 90th percentile OLS IV difference
absdiff = abs(IRF - IRFiv);
A = quantile(absdiff,.105);
B10 = absdiff <= A & absdiff >= A;

A = quantile(absdiff,.905);
B90 = absdiff <= A & absdiff >= A;


%%
f4=figure('Renderer', 'painters', 'Position', [10 10 700 300]);

A = polyval( polyfit(Ks,absdiff(B10),5), Ksi);
B1 = polyval( polyfit(Ks,2*IRFse(B10),5), Ksi);
B2 = polyval( polyfit(Ks,2*IRFivse(B10),5), Ksi);
B3 = polyval( polyfit(Ks,2*IRFdiffse(B10),5), Ksi);

%A = absdiff(B10);
%B1 = 2*IRFse(B10);
%B2 = 2*IRFivse(B10);
%B3 = 2*IRFdiffse(B10);

s1=subplot(2,2,1);
p1=plot(Ksi,A,"k",Ksi,B1,"--r",Ksi,B2,"--b",Ksi,B3,"--k")
ylim([.0 .25])
xticks([5 25 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)
title("10th percentile experiment")

A = polyval( polyfit(Ks,absdiff(B90),5), Ksi);
B1 = polyval( polyfit(Ks,2*IRFse(B90),5), Ksi);
B2 = polyval( polyfit(Ks,2*IRFivse(B90),5), Ksi);
B3 = polyval( polyfit(Ks,2*IRFdiffse(B90),5), Ksi);

s2=subplot(2,2,2);
p2=plot(Ksi,A,"k",Ksi,B1,"--r",Ksi,B2,"--b",Ksi,B3,"--k")
ylim([.0 .25])
xticks([5 25 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)
title("90th percentile experiment")

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'OLS$-$IV along estimated feedback direction','2 $\times$ OLS SE','2 $\times$ IV SE','2 $\times$ OLS$-$IV SE','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
axis(hLegend,'off');
box(leg,'off');
set(leg,'Position',posLegend);

% Adjusting positions of first two plots
poss1 = get(s1,'Position');
poss1(2) = poss1(2) - .37;
poss1(4) = poss1(4) + .37;
set(s1,'Position',poss1);

poss2 = get(s2,'Position');
poss2(2) = poss2(2) - .37;
poss2(4) = poss2(4) + .37;
set(s2,'Position',poss2);

%%

Demp = abs(COEF-COEFiv);
Temp_ols = abs(COEF) ./ SE;
Temp_iv = abs(COEFiv) ./ SEiv;
Temp = abs(COEF-COEFiv) ./ SE;
Temp_diff = abs(COEF-COEFiv) ./ SEDIFF;
Temp_diff2 = abs(COEF-COEFiv) ./ SE;
[max(Temp_ols,[],'all') max(Temp_iv,[],'all') max(Temp_diff,[],'all') max(Temp_diff2,[],'all')]

MDemp = zeros(1,10);
MTemp_diff2 = zeros(1,10);
MTemp_diff = zeros(1,10);
mTemp_ols = zeros(1,10);
mTemp_iv = zeros(1,10);
mTemp_diff = zeros(1,10);
for k = 1:Knum
    MDemp(k) = max(Demp(:,1:Ks(k),k),[],'all');
    MTemp_diff2(k) = max(Temp_diff2(:,1:Ks(k),k),[],'all');
    MTemp_diff(k) = max(Temp_diff(:,1:Ks(k),k),[],'all');
    mTemp_ols(k) = mean(Temp_ols(:,1:Ks(k),k) > 1.96,'all');
    mTemp_iv(k) = mean(Temp_iv(:,1:Ks(k),k) > 1.96,'all');
    mTemp_diff(k) = mean(Temp_diff(:,1:Ks(k),k) > 1.96,'all');
end
format bank
[MDemp; MTemp_diff2; MTemp_diff; mTemp_ols; mTemp_iv; mTemp_diff*100]


%%
% f5=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
% 
% A = quantile(COEF,.1);
% 
% s1=subplot(2,2,1);
% p1=plot(Ks,A)
% ylim([.0 4])
% xticks([5 25 50 100])
% xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)
% title("10th percentile experiment")
% 
% A = quantile(COEF,.9);
% 
% s2=subplot(2,2,2);
% p1=plot(Ks,A)
% ylim([.0 4])
% xticks([5 25 50 100])
% xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)
% title("90th percentile experiment")
% 
% % Plotting the legend
% hLegend = subplot(2,2,3.5);
% posLegend = get(hLegend,'Position');
% posLegend(2) = posLegend(2) - .23; 
% leg=legend(hLegend,p1,'Maximum coefficient in abs(OLS-IV)/se(OLS-IV) ','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
% axis(hLegend,'off');
% box(leg,'off');
% set(leg,'Position',posLegend);
% 
% % Adjusting positions of first two plots
% poss1 = get(s1,'Position');
% poss1(2) = poss1(2) - .37;
% poss1(4) = poss1(4) + .37;
% set(s1,'Position',poss1);
% 
% poss2 = get(s2,'Position');
% poss2(2) = poss2(2) - .37;
% poss2(4) = poss2(4) + .37;
% set(s2,'Position',poss2);


%%
print(f1,'Figures\BiasFigure_Fred','-depsc','-tiff')
print(f2,'Figures\SizeFigure_Fred','-depsc','-tiff')
print(f4,'Figures\DiffFigureFeedback_Fred','-depsc','-tiff')
%print(f5,'Figures\DiffFigureMax_Fred','-depsc','-tiff')