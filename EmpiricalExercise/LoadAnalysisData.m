%% Load data
function [analysisData,analysisData_names] = LoadAnalysisData

small = 1.0e-10;
big = 1.0e+6;

% ----------- Sample Period, Calendars and so forth
[dnobs_m,calvec_m,calds_m] = calendar_make([1959 1],[2014 12],12);  % Monthly Calendar
[dnobs_q,calvec_q,calds_q] = calendar_make([1959 1],[2014 4],4);    % Quarterly Calendar

% -- Load Data
load_data=1;
% Demeaning Parameters
i_demean = 1;  % 0 Do Nothing
% 1 Eliminate low-frequency by local Demeaning
% 2 Eliminate low-frequency trend by full-sample demeaning;

% Transformations
take_differences = 2; 
% 1 takes first and second differences of nearly non-stationary sequences
% 2 uses Hamilton (2018) to transform the data
% 0 does nothing


bw_bw = 100;   % Bi-Weight Parameter for local demeaning
datain_all;      % datain_all reads in the full dataset .. all variables, etc. saved in datain.xx

% Sampling parameters
smpl_par.nfirst = [1963 1];       % start date
smpl_par.nlast  = [2013 4];       % end date
smpl_par.calvec = datain.calvec;  % calendar
smpl_par.nper   = 4;              % number of periods a year

%% Estimation data
size(datain.bpdata); % Contains 224 observations and 207 variables
est_data = datain.bpdata(:,datain.bpinclcode==1); % Contains 224 observations and 139 variables
est_names = datain.bplabvec_long(:,datain.bpinclcode==1);


%% Sample period
[istart, iend] = smpl_HO(smpl_par);
istart = max(istart,1);
iend = min(iend,size(est_data,1));

% Extract data from sampling period
xdata = est_data(istart:iend,:); % Contains 204 observations and 111 variables
xdataSize = size(xdata);
miss = sum(isnan(xdata),1)/xdataSize(1); %Fraction of missings per column

% Drop variables with missing values
analysisData = xdata(:,miss==0); % Contains 204 observations and 111 variables
analysisData_names = est_names(:,miss==0)';

end