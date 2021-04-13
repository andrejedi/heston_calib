S0=100;
K=140;
T=1;
q=0;
alpha=1.51;
r=0.1;
v0=0.5^2;
kappa=1;
theta=0.5^2;
sigma=0.5;
rho=0.5;

N=5;

 CallPrice_v1  = tylkozadzialaj( S0,K,T,q,alpha,r,v0,kappa,theta,sigma,rho,N );

I1=integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),0,N)+...
    integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),N,5*N)+...
    integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),5*N,10*N);
I2=integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-N,0)+...
    integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-5*N,-N)+...
    integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-10*N,-5*N);
CallPrice_v1=S0*exp(-r*T)/(2*pi)*real(I1+I2)

CallPrice_v2=-S0*exp(-r*T)/(2*pi)*(integral(@(w) (K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)+...
    (K/S0).^(1-alpha-1i*w)./((w-1i*(alpha-1)).*(w-1i*alpha)).*FT_prob_dens(w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),0,N)+...
    integral(@(w) (K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)+...
    (K/S0).^(1-alpha-1i*w)./((w-1i*(alpha-1)).*(w-1i*alpha)).*FT_prob_dens(w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),N,5*N)+...
    integral(@(w) (K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)+...
    (K/S0).^(1-alpha-1i*w)./((w-1i*(alpha-1)).*(w-1i*alpha)).*FT_prob_dens(w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),5*N,10*N))

A=log(S0/K)+(r-q)*T-rho*(theta*kappa*T+v0)/sigma;
B=(theta*kappa*T+v0)/sigma*sqrt(1-rho^2);
phi_star1=atan(-A/B);
phi_star=max(phi_star1,-pi/3);
phi_star=phi_star1;
ratio=sqrt(A^2+B^2)/B;
acc1=10^(-12);
acc2=10^(-14);
N=1.31;
II1=integral(@(R) (K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star),0,N,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star),N,5*N,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star),5*N,10*N,'RelTol',acc1,'AbsTol',acc2);
II2=integral(@(R) (K/S0).^(1-alpha-1i*R*exp(-1i*phi_star))./((-R*exp(-1i*phi_star)+1i*(alpha-1)).*(-R*exp(-1i*phi_star)+1i*alpha)).*FT_prob_dens(R*exp(-1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(-1i*phi_star),0,N,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha-1i*R*exp(-1i*phi_star))./((-R*exp(-1i*phi_star)+1i*(alpha-1)).*(-R*exp(-1i*phi_star)+1i*alpha)).*FT_prob_dens(R*exp(-1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(-1i*phi_star),N,5*N,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha-1i*R*exp(-1i*phi_star))./((-R*exp(-1i*phi_star)+1i*(alpha-1)).*(-R*exp(-1i*phi_star)+1i*alpha)).*FT_prob_dens(R*exp(-1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(-1i*phi_star),5*N,10*N,'RelTol',acc1,'AbsTol',acc2);
CallPrice_v3=-S0*exp(-r*T)/(2*pi)*real(II1+II2)

IV1=blsimpv(S0,K,r,T,CallPrice_v1)
IV2=blsimpv(S0,K,r,T,CallPrice_v2)
IV3=blsimpv(S0,K,r,T,CallPrice_v3)

R=0.001:0.001:50;
% phi_star=0;
integrand1=(K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star);
zz=imag(integrand1)./abs(integrand1);
plot(R,asin(zz))
phi_star=-0+0.0;
integrand2=(K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star);
plot(R,[abs(integrand1);abs(integrand2)])
legend('phi star','-pi/4')
zzz=abs(integrand1)-abs(integrand2);
