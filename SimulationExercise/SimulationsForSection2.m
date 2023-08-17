%% Setup
clear

% Global Parameters
T = 200;
rep = 2000;

% Weight on feedback mechanism
a = 1.5;
% Number of feedback periods in the data
Sdata = 1;
% Number of leads used in the IV
S = 1;
% Number of periods to assess potential for bias with
Sb = 20;
% IV type: 'lead' or 'future' (future use entire future of regressors for the IV)
type = 'lead';
% DGP for strictly exongeous regressors: 'AR' or 'MA'
dgp = 'AR';

% Where to save and whether or not to save (saving=1 is saving)
saving = 1;
f1Name = 'Figures\biasFigure';
f2Name = 'Figures\biasFigureOLSIV';
f3Name = 'Figures\biasFigureOLS';
f31Name = 'Figures\biasFigureOLSTr';
f4Name = 'Figures\biasBoundsFigure';
f5Name = 'Figures\biasBoundsFigureLargeL';
f6Name = 'Figures\biasBoundsFigureOLS';
f7Name = 'Figures\biasBoundsFigureIV';
workspaceName = 'SimulationsForSection2.mat';

%% Simulations varying number of regressors (K)

% Setup for this simulation
% Number of parameters
Kvec = [4,6,8,10,12,14,16,18,20,30,40,50,60,70,80,90,100];
lK = length(Kvec);
% Autocorrelation coefficient for strictly exogenous regressors
rho = 0.8;

tic
summaryK = cell(lK,1);
resK = cell(lK,1);
parfor k = 1:lK
    % Fix seed for random number generation
    rng(k);

    K = Kvec(k);
    resK{k} = zeros(rep,5+2*Sb);
    for i=1:rep
        % Generate data
        [y,X] = genData(T,K,rho,a,Sdata,dgp);
        % First and second coefficients from OLS and IV
        [ols,iv,tr,Bbols,Bbiv] = reg_func(y,X,S,type,1:2,Sb);
        resK{k}(i,:) = [ols(1),ols(2),iv(1),iv(2),tr{1},Bbols',Bbiv'];
    end
    % Compute summary statistics from the simulation
    summaryK{k} = [mean(resK{k});sqrt( var(resK{k}) )];
end
toc

%% Simulations varying the autocorrelation in regressors (rho)

% Setup for this simulation
% Number of parameters
K = 50;
% Autocorrelation coefficient for strictly exogenous regressors
RHOvec = [0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.92,0.94,0.96,0.98];
lRHO = length(RHOvec);


tic
summaryRHO = cell(lRHO,1);
resRHO = cell(lRHO,1);
parfor k = 1:lRHO
    % Fix seed for random number generation
    rng(100+k);

    rho = RHOvec(k);
    resRHO{k} = zeros(rep,5+2*Sb);
    for i=1:rep
        % Generate data
        [y,X] = genData(T,K,rho,a,Sdata,dgp);
        % First and second coefficients from OLS and IV
        [ols,iv,tr,Bbols,Bbiv] = reg_func(y,X,S,type,1:2,Sb);
        resRHO{k}(i,:) = [ols(1),ols(2),iv(1),iv(2),tr{1},Bbols',Bbiv'];
    end
    % Compute summary statistics from the simulation
    summaryRHO{k} = [ mean(resRHO{k}); sqrt( var(resRHO{k}) )];
end
toc

%% Plotting the results

%% Preparations
% Varying K (first coordinate of hat beta)
SumK = zeros(lK,5+2*Sb);
for k=1:lK
    SumK(k,1) = summaryK{k}(1,1) ;
    SumK(k,2)  = summaryK{k}(2,1) ;
    SumK(k,3)  = summaryK{k}(1,3) ;
    SumK(k,4)   = summaryK{k}(2,3) ;
    SumK(k,5)   = summaryK{k}(1,5) / T ;
    for j=1:Sb
        SumK(k,5+j) = summaryK{k}(1,5+j);
        SumK(k,5+Sb+j)  = summaryK{k}(1,5+Sb+j);
    end
end

Kveci = Kvec(1):0.1:Kvec(lK);
SumKi = zeros(5+2*Sb,length(Kveci));
for k = 1:5+2*Sb
    SumKi(k,:) = polyval( polyfit(Kvec,SumK(:,k),6), Kveci);
end

% Varying rho
SumRHO = zeros(lRHO,5+2*Sb);
for k=1:lRHO
    SumRHO(k,1) = summaryRHO{k}(1,1) ;
    SumRHO(k,2)  = summaryRHO{k}(2,1) ;
    SumRHO(k,3)  = summaryRHO{k}(1,3) ;
    SumRHO(k,4)   = summaryRHO{k}(2,3) ;
    SumRHO(k,5)   = summaryRHO{k}(1,5) / T ;
    for j=1:Sb
        SumRHO(k,5+j) = summaryRHO{k}(1,5+j);
        SumRHO(k,5+Sb+j)  = summaryRHO{k}(1,5+Sb+j);
    end
end

RHOveci = RHOvec(1):0.001:RHOvec(lRHO);
SumRHOi = zeros(5+2*Sb,length(RHOveci));
for k = 1:5+2*Sb
    SumRHOi(k,:) = polyval( polyfit(RHOvec,SumRHO(:,k),6), RHOveci);
end


%% Plotting all results
f1=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
    
% Plotting for K
s1=subplot(2,2,1);
p1=plot(Kveci,abs(SumKi(1,:)),'r',Kveci,SumKi(2,:),'--r',Kveci,abs(SumKi(3,:)),'b',Kveci,SumKi(4,:),'--b',Kveci,abs(SumKi(5,:)),':k');
ylim([-.01 .17])
xticks([4 20 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)

% Plotting for rho
s2=subplot(2,2,2);
p2=plot(RHOveci,abs(SumRHOi(1,:)),'r',RHOveci,SumRHOi(2,:),'--r',RHOveci,abs(SumRHOi(3,:)),'b',RHOveci,SumRHOi(4,:),'--b',RHOveci,abs(SumRHOi(5,:)),':k');
ylim([-.01 .17])
xticks([0 .3 .5 .98])
xlabel('$\rho$: Autocorrelation Coefficient','Interpreter','latex', 'FontSize', 12)

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


%% Plotting OLS and IV
f2=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
    
% Plotting for K
s1=subplot(2,2,1);
p1=plot(Kveci,abs(SumKi(1,:)),'r',Kveci,SumKi(2,:),'--r',Kveci,abs(SumKi(3,:)),'b',Kveci,SumKi(4,:),'--b');
ylim([-.01 .17])
xticks([4 20 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)

% Plotting for rho
s2=subplot(2,2,2);
p2=plot(RHOveci,abs(SumRHOi(1,:)),'r',RHOveci,SumRHOi(2,:),'--r',RHOveci,abs(SumRHOi(3,:)),'b',RHOveci,SumRHOi(4,:),'--b');
ylim([-.01 .17])
xticks([0 .3 .5 .98])
xlabel('$\rho$: Autocorrelation Coefficient','Interpreter','latex', 'FontSize', 12)

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Bias (OLS)','Stddev (OLS)','Bias (IV)','Stddev (IV)','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
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

%% Plotting OLS only
f3=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
    
% Plotting for K
s1=subplot(2,2,1);
p1=plot(Kveci,abs(SumKi(1,:)),'r',Kveci,SumKi(2,:),'--r');
ylim([-.01 .17])
xticks([4 20 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)

% Plotting for rho
s2=subplot(2,2,2);
p2=plot(RHOveci,abs(SumRHOi(1,:)),'r',RHOveci,SumRHOi(2,:),'--r');
ylim([-.01 .17])
xticks([0 .3 .5 .98])
xlabel('$\rho$: Autocorrelation Coefficient','Interpreter','latex', 'FontSize', 12)

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Bias (OLS)','Stddev (OLS)','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
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


%% Plotting OLS with trace
f31=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
    
% Plotting for K
s1=subplot(2,2,1);
p1=plot(Kveci,abs(SumKi(1,:)),'r',Kveci,SumKi(2,:),'--r',Kveci,abs(SumKi(5,:)),':k');
ylim([-.01 .17])
xticks([4 20 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)

% Plotting for rho
s2=subplot(2,2,2);
p2=plot(RHOveci,abs(SumRHOi(1,:)),'r',RHOveci,SumRHOi(2,:),'--r',RHOveci,abs(SumRHOi(5,:)),':k');
ylim([-.01 .17])
xticks([0 .3 .5 .98])
xlabel('$\rho$: Autocorrelation Coefficient','Interpreter','latex', 'FontSize', 12)

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Bias (OLS)','Stddev (OLS)','Trace/$T$','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
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


%% Plotting Bias bounds
f4=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
    
% Plotting for K
s1=subplot(2,2,1);
p1=plot(Kveci,SumKi(5+1,:).*SumKi(2,:),'r',Kveci,SumKi(5+6,:).*SumKi(2,:),':r',Kveci,SumKi(5+20,:).*SumKi(2,:),'-.r',Kveci,SumKi(5+Sb+1,:).*SumKi(4,:),'b',Kveci,SumKi(5+Sb+6,:).*SumKi(4,:),':b',Kveci,SumKi(5+Sb+20,:).*SumKi(4,:),'-.b');
ylim([-.01 .25])
xticks([4 20 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)

% Plotting for rho
s2=subplot(2,2,2);
p2=plot(RHOveci,SumRHOi(5+1,:).*SumRHOi(2,:),'r',RHOveci,SumRHOi(5+6,:).*SumRHOi(2,:),':r',RHOveci,SumRHOi(5+20,:).*SumRHOi(2,:),'-.r',RHOveci,SumRHOi(5+Sb+1,:).*SumRHOi(4,:),'b',RHOveci,SumRHOi(5+Sb+6,:).*SumRHOi(4,:),':b',RHOveci,SumRHOi(5+Sb+20,:).*SumRHOi(4,:),'-.b');
ylim([-.01 .25])
xticks([0 .3 .5 .98])
xlabel('$\rho$: Autocorrelation Coefficient','Interpreter','latex', 'FontSize', 12)

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Bound$_1$ (OLS)','Bound$_6$ (OLS)','Bound$_{20}$ (OLS)','Bound$_1$ (IV)','Bound$_6$ (IV)','Bound$_{20}$ (IV)','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
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


%% Plotting Bias bounds for large L only
f5=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
    
% Plotting for K
s1=subplot(2,2,1);
p1=plot(Kveci,SumKi(5+6,:).*SumKi(2,:),':r',Kveci,SumKi(5+20,:).*SumKi(2,:),'-.r',Kveci,SumKi(5+Sb+6,:).*SumKi(4,:),':b',Kveci,SumKi(5+Sb+20,:).*SumKi(4,:),'-.b');
ylim([-.01 .25])
xticks([4 20 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)

% Plotting for rho
s2=subplot(2,2,2);
p2=plot(RHOveci,SumRHOi(5+6,:).*SumRHOi(2,:),'r',RHOveci,SumRHOi(5+20,:).*SumRHOi(2,:),'-.r',RHOveci,SumRHOi(5+Sb+6,:).*SumRHOi(4,:),':b',RHOveci,SumRHOi(5+Sb+20,:).*SumRHOi(4,:),'-.b');
ylim([-.01 .25])
xticks([0 .3 .5 .98])
xlabel('$\rho$: Autocorrelation Coefficient','Interpreter','latex', 'FontSize', 12)

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Bound$_6$ (OLS)','Bound$_{20}$ (OLS)','Bound$_6$ (IV)','Bound$_{20}$ (IV)','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
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


%% Plotting Bias bounds for OLS only
f6=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
    
% Plotting for K
s1=subplot(2,2,1);
p1=plot(Kveci,SumKi(5+1,:).*SumKi(2,:),'r',Kveci,SumKi(5+2,:).*SumKi(2,:),'--r',Kveci,SumKi(5+4,:).*SumKi(2,:),'--r',Kveci,SumKi(5+6,:).*SumKi(2,:),':r',Kveci,SumKi(5+20,:).*SumKi(2,:),'-.r');
ylim([-.01 .25])
xticks([4 20 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)

% Plotting for rho
s2=subplot(2,2,2);
p2=plot(RHOveci,SumRHOi(5+1,:).*SumRHOi(2,:),'r',RHOveci,SumRHOi(5+2,:).*SumRHOi(2,:),'--r',RHOveci,SumRHOi(5+4,:).*SumRHOi(2,:),'--r',RHOveci,SumRHOi(5+6,:).*SumRHOi(2,:),':r',RHOveci,SumRHOi(5+20,:).*SumRHOi(2,:),'-.r');
ylim([-.01 .25])
xticks([0 .3 .5 .98])
xlabel('$\rho$: Autocorrelation Coefficient','Interpreter','latex', 'FontSize', 12)

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Bound$_1$ (OLS)','Bound$_{2}$ (OLS)','Bound$_4$ (OLS)','Bound$_6$ (OLS)','Bound$_{20}$ (OLS)','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
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

%% Plotting Bias bounds for IV only
f7=figure('Renderer', 'painters', 'Position', [10 10 700 300]);
    
% Plotting for K
s1=subplot(2,2,1);
p1=plot(Kveci,SumKi(5+1+Sb,:).*SumKi(4,:),'b',Kveci,SumKi(5+2+Sb,:).*SumKi(4,:),'--b',Kveci,SumKi(5+4+Sb,:).*SumKi(4,:),'--b',Kveci,SumKi(5+6+Sb,:).*SumKi(4,:),':b',Kveci,SumKi(5+20+Sb,:).*SumKi(4,:),'-.b');
ylim([-.01 .25])
xticks([4 20 50 100])
xlabel('$K$: Number of Regressors','Interpreter','latex', 'FontSize', 12)

% Plotting for rho
s2=subplot(2,2,2);
p2=plot(RHOveci,SumRHOi(5+1+Sb,:).*SumRHOi(4,:),'b',RHOveci,SumRHOi(5+2+Sb,:).*SumRHOi(4,:),'--b',RHOveci,SumRHOi(5+4+Sb,:).*SumRHOi(4,:),'--b',RHOveci,SumRHOi(5+6+Sb,:).*SumRHOi(4,:),':b',RHOveci,SumRHOi(5+20+Sb,:).*SumRHOi(4,:),'-.b');
ylim([-.01 .25])
xticks([0 .3 .5 .98])
xlabel('$\rho$: Autocorrelation Coefficient','Interpreter','latex', 'FontSize', 12)

% Plotting the legend
hLegend = subplot(2,2,3.5);
posLegend = get(hLegend,'Position');
posLegend(2) = posLegend(2) - .23; 
leg=legend(hLegend,p1,'Bound$_1$ (OLS)','Bound$_{2}$ (OLS)','Bound$_4$ (OLS)','Bound$_6$ (OLS)','Bound$_{20}$ (OLS)','Interpreter','latex','Orientation','horizontal', 'FontSize', 12);
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

%% Save Figures and Workspace to file
if saving == 1
    % Save Figures to file
    print(f1,f1Name,'-depsc','-tiff')
    print(f2,f2Name,'-depsc','-tiff')
    print(f3,f3Name,'-depsc','-tiff')
    print(f31,f31Name,'-depsc','-tiff')
    print(f4,f4Name,'-depsc','-tiff')
    print(f5,f5Name,'-depsc','-tiff')
    print(f6,f6Name,'-depsc','-tiff')
    print(f7,f7Name,'-depsc','-tiff')

    % Save Workspace to file
    close all
    save(fullfile('SavedResults',workspaceName));
end
