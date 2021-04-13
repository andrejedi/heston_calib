%% ONE-BY-ONE calibration 

underlyings = ['ALKS' 'AMAG' 'VOC' 'WNC' 'ZN'];

CallPrices = readtable('ZN.csv');
MATLABDate = x2mdate(CallPrices.Expiration,0,'datetime');
CallPrices.Expiration = MATLABDate;
t0 = datenum('01/01/2015');
CallPrices.TimeToMaturity = (datenum(CallPrices.Expiration) - t0)/360;

S0=CallPrices.UnderlyingPrice;
MarketPrice=CallPrices.Ask;
K=CallPrices.Strike;
T=CallPrices.TimeToMaturity;

q=0;
alpha=1.51;
r=0.1;
N=100;

options = optimset('MaxFunEvals',50000);
params = array2table(zeros(size(K,1),5));
startparameters = [0.02 1 0 0.021 1];

for i=1:numel(K)
    MP = MarketPrice(i);
    S0i = S0(i);
    Ki = K(i);
    Ti = T(i);
    MarketIVi = MarketIV(i);
    calibrator=@(x)(MP-heston_call(S0i,Ki,Ti,q,alpha,r,x(1),x(2),x(3),x(4),x(5),N))^2/MarketIVi^2;
    params(i,:) = array2table(fminsearch(calibrator,startparameters,options));
end

%% MULTI NON LINEAR OPTIMALIZATION

underlyings = ['ALKS' 'AMAG' 'VOC' 'WNC' 'ZN'];

CallPrices = readtable('ALKS.csv');
MATLABDate = x2mdate(CallPrices.Expiration,0,'datetime');
CallPrices.Expiration = MATLABDate;
t0 = datenum('01/01/2015');
CallPrices.TimeToMaturity = (datenum(CallPrices.Expiration) - t0)/360;

S0=CallPrices.UnderlyingPrice;
MarketPrice=CallPrices.Ask;
K=CallPrices.Strike;
T=CallPrices.TimeToMaturity;

q=0;
alpha=1.51;
r=0.1;
N=100;

MarketIV = impliedvol_c(MarketPrice,S0,K,T,r);

calibrator=@(x)(MarketPrice-heston_call_c(S0,K,T,q,alpha,r,x(1),x(2),x(3),x(4),x(5),N)).^2./MarketIV.^2;
startparameters = [0.01 0.01 0.01 0.01 0.01];
xopt = lsqnonlin(calibrator,startparameters);

%,[0 0 0 0 -1],[1 100 1 1 1]

%% PRICE VALIDATION

underlyings = ['ALKS' 'AMAG' 'VOC' 'WNC' 'ZN'];

CallPrices = readtable('ZN.csv');
MATLABDate = x2mdate(CallPrices.Expiration,0,'datetime');
CallPrices.Expiration = MATLABDate;
t0 = datenum('01/01/2015');
CallPrices.TimeToMaturity = (datenum(CallPrices.Expiration) - t0)/360;

S0=CallPrices.UnderlyingPrice;
K=CallPrices.Strike;
T=CallPrices.TimeToMaturity;

q=0;
alpha=1.51;
r=0.1;
N=100;

v0=	0.35726745	;
kappa=	0.70918595	;
theta=	0.316856688	;
sigma=	0.933644321	;
rho=	0.171654455	;

CallPrice_v1_c = heston_call_c(S0,K,T,q,alpha,r,v0,kappa,theta,sigma,rho,N);

%% alpha calibration - on-eby-one
underlyings = ['ALKS' 'AMAG' 'VOC' 'WNC' 'ZN'];

CallPrices = readtable('ZN.csv');
MATLABDate = x2mdate(CallPrices.Expiration,0,'datetime');
CallPrices.Expiration = MATLABDate;
t0 = datenum('01/01/2015');
CallPrices.TimeToMaturity = (datenum(CallPrices.Expiration) - t0)/360;

S0=CallPrices.UnderlyingPrice;
MarketPrice=CallPrices.Ask;
K=CallPrices.Strike;
T=CallPrices.TimeToMaturity;
N=100;

MarketIV = impliedvol_c(MarketPrice,S0,K,T,r);
options = optimset('MaxFunEvals',50000);
params = array2table(zeros(size(K,1),6));
startparameters = [0, 0.01 0.01 0 0.01 0.01];

for i=1:numel(K)
    MP = MarketPrice(i);
    S0i = S0(i);
    Ki = K(i);
    Ti = T(i);
    MarketIVi = MarketIV(i);
    calibrator=@(x)(MP-heston_call(S0i,Ki,Ti,q,x(1),r,x(2),x(3),x(4),x(5),x(6),N))^2/MarketIVi^2;
    params(i,:) = array2table(fminsearch(calibrator,startparameters,options));
end

%% MULTI NON LINEAR OPTIMALIZATION

underlyings = ['ALKS' 'AMAG' 'VOC' 'WNC' 'ZN'];

CallPrices = readtable('ALKS.csv');
MATLABDate = x2mdate(CallPrices.Expiration,0,'datetime');
CallPrices.Expiration = MATLABDate;
t0 = datenum('01/01/2015');
CallPrices.TimeToMaturity = (datenum(CallPrices.Expiration) - t0)/360;

S0=CallPrices.UnderlyingPrice;
MarketPrice=CallPrices.Ask;
K=CallPrices.Strike;
T=CallPrices.TimeToMaturity;

q=0;
r=0.1;
N=100;

MarketIV = impliedvol_c(MarketPrice,S0,K,T,r);

calibrator=@(x)(MarketPrice-heston_call_c(S0,K,T,q,x(1),r,x(2),x(3),x(4),x(5),x(6),N)).^2./MarketIV.^2;
startparameters = [0.01 0.01 0.01 0.01 0.01 0.01];
xopt = lsqnonlin(calibrator,startparameters,[0 0 0 0 0 -1],[500 1 100 1 1 1]);

%,[0 0 0 0 0 -1],[500 1 100 1 1 1]

%% PRICE VALIDATION (with calibrated alpha)

underlyings = ['ALKS' 'AMAG' 'VOC' 'WNC' 'ZN'];

CallPrices = readtable('ZN.csv');
MATLABDate = x2mdate(CallPrices.Expiration,0,'datetime');
CallPrices.Expiration = MATLABDate;
t0 = datenum('01/01/2015');
CallPrices.TimeToMaturity = (datenum(CallPrices.Expiration) - t0)/360;

S0=CallPrices.UnderlyingPrice;
K=CallPrices.Strike;
T=CallPrices.TimeToMaturity;

q=0;
r=0.1;
N=100;
alpha=	161.2462102	;
v0=	0.902823167	;
kappa=	0.696194239	;
theta=	0.893247286	;
sigma=	0.858985501	;
rho=	0.905850343	;


CallPrice_v1_c = heston_call_c(S0,K,T,q,alpha,r,v0,kappa,theta,sigma,rho,N);