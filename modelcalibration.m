
%% READ MARKET DATA
CallPrices = readtable('CallPrice.csv');

MATLABDate = x2mdate(CallPrices.Expiration,0,'datetime');
CallPrices.Expiration = MATLABDate;

t0 = datenum('01/01/2018');
CallPrices.TimeToMaturity = (datenum(CallPrices.Expiration) - t0)/360;

%% ONE OPTION CALIBRATION EXAMPLE
S0=CallPrices.UnderlyingPrice(1);
K=CallPrices.Strike(1);
T=CallPrices.TimeToMaturity(1);
MarketPrice=CallPrices.Ask(1);

q=0;
alpha=1.51;
r=0.1;
v0=0.5^2;
kappa=1;
theta=0.5^2;
sigma=0.5;
rho=0.5;
N=5;

CallPrice_v1 = heston_call( S0,K,T,q,alpha,r,v0,kappa,theta,sigma,rho,N);
MarketIV = blsimpv(S0,K,r,T,MarketPrice);

options = optimset('MaxFunEvals',5000);
startparameters = [0.02 1 -0 0.021 1];
calibrator=@(x)(MarketPrice-heston_call(S0,K,T,q,alpha,r,x(1),x(2),x(3),x(4),x(5),N))^2/MarketIV^2;
fminsearch(calibrator,startparameters, options)
 % output = [0.0297    0.5681   -0.0204    0.0093    4.7924]


%% MULTI NON LINEAR OPTIMALIZATION

S0=CallPrices.UnderlyingPrice;
MarketPrice=CallPrices.Ask;
K=CallPrices.Strike;
T=CallPrices.TimeToMaturity;

q=0;
alpha=1.51;
r=0.1;
v0=0.5^2;
kappa=1;
theta=0.5^2;
sigma=0.5;
rho=0.5;
N=5;


CallPrice_v1_c = heston_call_c( S0,K,T,q,alpha,r,v0,kappa,theta,sigma,rho,N );
MarketIV = impliedvol_c(MarketPrice,S0,K,T,r);

calibrator=@(x)(MarketPrice-heston_call_c(S0,K,T,q,alpha,r,x(1),x(2),x(3),x(4),x(5),N)).^2./MarketIV.^2;
startparameters = [0.0297    0.5681   -0.0204    0.0093    4.7924];

options = optimoptions('lsqnonlin', 'Display', 'iter');
xopt = lsqnonlin(calibrator,startparameters);


%% LOOP CHECK
options = optimset('MaxFunEvals',5000);
params = array2table(zeros(size(K,1),5));
startparameters = [0.02 1 -0 0.021 1];

calibrator=@(x)(MP-heston_call(S0i,Ki,Ti,q,alpha,r,x(1),x(2),x(3),x(4),x(5),N))^2/MarketIVi^2;
for i=1:numel(K)
    MP = MarketPrice(i);
    S0i = S0(i);
    Ki = K(i);
    Ti = T(i);
    MarketIVi = MarketIV(i);
    calibrator=@(x)(MP-heston_call(S0i,Ki,Ti,q,alpha,r,x(1),x(2),x(3),x(4),x(5),N))^2/MarketIVi^2;
    params(i,:) = array2table(fminsearch(calibrator,startparameters,options));
end