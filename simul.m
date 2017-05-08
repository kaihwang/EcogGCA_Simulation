%% simulate timeseries, for real analysis this part will be reading in data

nvar = 3; %3 nodes in the network, #1:frontal, #2:thalamus, #3:parietal
p = 20; %maximum model order, assuming 50hz sampling
amo = 5; %acutual time lag order where interactions happen, here assuming it is alpha, 1000/50 = 5
fs        = 50;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

%setup fake var model A, assuming interactions happening at alpha frequency, which
%would be 5
A = zeros(nvar,nvar,p);
A(2,1,4:6) = 1; %from frontal to thalamus at 4-6 lags (alpha freq)
A(3,2,4:6) = 1; %from thalamus to pareital at 4-6 lags (alpha freq)

%residual cov matrix
SIGT = eye(nvar);

%simulate TS
nobs = 1000; %number of observation(timepoints), here assume 1000
ntrials = 100; %number of trials, here assume 100
X = var_to_tsdata(A,SIGT,nobs,ntrials); %X is the data variable, nvar(ROIs) x time x trials

%% estimate model order

[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,p);


figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
fprintf('actual model order     = %d\n',amo);
morder = moBIC; %use optimal order according to BIC

%% do GCA 
[AR,SIG] = tsdata_to_var(X,morder); %estimate VAR model
[G,info] = var_to_autocov(AR,SIG,p); %autocovaraince 
f = autocov_to_spwcgc(G,fres);

figure(3); clf;
plot_spw(f,fs);
