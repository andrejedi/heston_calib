clear all
clc
S0=100;
T=0.05;
q=0;
alpha=1.5;
r=0.1;
v0=0.5^2;
kappa=0.1;
theta=0.9^2;
sigma=0.5;
rho=-0.5;

acc1=10^(-13);
acc2=10^(-15);

for j=1:2100
    N=100;
    KK(j)=j*0.1;
    K=KK(j);
    I1=integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),0,N/10,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),N/10,N,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),N,5*N,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),5*N,10*N,'RelTol',acc1,'AbsTol',acc2);
    I2=integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-N/10,0,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-N,-N/10,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-5*N,-N,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-10*N,-5*N,'RelTol',acc1,'AbsTol',acc2);
    CallPrice_v1(j)=S0*exp(-r*T)/(2*pi)*real(I1+I2);

    CallPrice_v2(j)=-S0*exp(-r*T)/(2*pi)*(integral(@(w) (K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)+...
        (K/S0).^(1-alpha-1i*w)./((w-1i*(alpha-1)).*(w-1i*alpha)).*FT_prob_dens(w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),0,N/10,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) (K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)+...
        (K/S0).^(1-alpha-1i*w)./((w-1i*(alpha-1)).*(w-1i*alpha)).*FT_prob_dens(w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),N/10,N,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) (K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)+...
        (K/S0).^(1-alpha-1i*w)./((w-1i*(alpha-1)).*(w-1i*alpha)).*FT_prob_dens(w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),N,5*N,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) (K/S0).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)+...
        (K/S0).^(1-alpha-1i*w)./((w-1i*(alpha-1)).*(w-1i*alpha)).*FT_prob_dens(w-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),5*N,10*N,'RelTol',acc1,'AbsTol',acc2));
    
    alpha1=0.1;
    I1=integral(@(w) -(K/S0).^(1-alpha+1i*(w-(alpha+alpha1)*1i))./(((w-(alpha+alpha1)*1i)+1i*(alpha-1)).*((w-(alpha+alpha1)*1i)+1i*alpha)).*FT_prob_dens(-(w-(alpha+alpha1)*1i)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),0,N/10,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*(w-(alpha+alpha1)*1i))./(((w-(alpha+alpha1)*1i)+1i*(alpha-1)).*((w-(alpha+alpha1)*1i)+1i*alpha)).*FT_prob_dens(-(w-(alpha+alpha1)*1i)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),N/10,N,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*(w-(alpha+alpha1)*1i))./(((w-(alpha+alpha1)*1i)+1i*(alpha-1)).*((w-(alpha+alpha1)*1i)+1i*alpha)).*FT_prob_dens(-(w-(alpha+alpha1)*1i)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),N,5*N,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*(w-(alpha+alpha1)*1i))./(((w-(alpha+alpha1)*1i)+1i*(alpha-1)).*((w-(alpha+alpha1)*1i)+1i*alpha)).*FT_prob_dens(-(w-(alpha+alpha1)*1i)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),5*N,10*N,'RelTol',acc1,'AbsTol',acc2);
    I2=integral(@(w) -(K/S0).^(1-alpha+1i*(w-(alpha+alpha1)*1i))./(((w-(alpha+alpha1)*1i)+1i*(alpha-1)).*((w-(alpha+alpha1)*1i)+1i*alpha)).*FT_prob_dens(-(w-(alpha+alpha1)*1i)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-N/10,0,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*(w-(alpha+alpha1)*1i))./(((w-(alpha+alpha1)*1i)+1i*(alpha-1)).*((w-(alpha+alpha1)*1i)+1i*alpha)).*FT_prob_dens(-(w-(alpha+alpha1)*1i)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-N,-N/10,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*(w-(alpha+alpha1)*1i))./(((w-(alpha+alpha1)*1i)+1i*(alpha-1)).*((w-(alpha+alpha1)*1i)+1i*alpha)).*FT_prob_dens(-(w-(alpha+alpha1)*1i)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-5*N,-N,'RelTol',acc1,'AbsTol',acc2)+...
        integral(@(w) -(K/S0).^(1-alpha+1i*(w-(alpha+alpha1)*1i))./(((w-(alpha+alpha1)*1i)+1i*(alpha-1)).*((w-(alpha+alpha1)*1i)+1i*alpha)).*FT_prob_dens(-(w-(alpha+alpha1)*1i)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho),-10*N,-5*N,'RelTol',acc1,'AbsTol',acc2);
    CallPrice_v4(j)=S0-K*exp(-r*T)+S0*exp(-r*T)/(2*pi)*real(I1+I2);
    
    
    A=log(S0/K)+(r-q)*T-rho*(theta*kappa*T+v0)/sigma;
    B=(theta*kappa*T+v0)/sigma*sqrt(1-rho^2);
    phi_star1(j)=atan(-A/B);
    phi_star=max(phi_star1(j),-pi/4);
%     phi_star=phi_star1(j);
    ratio=sqrt(A^2+B^2)/B;
    N=2;
    II1=integral(@(R) (K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star),0,N/3,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star),N/3, N,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star),N,5*N,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha+1i*R*exp(1i*phi_star))./((R*exp(1i*phi_star)+1i*(alpha-1)).*(R*exp(1i*phi_star)+1i*alpha)).*FT_prob_dens(-R*exp(1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(1i*phi_star),5*N,10*N,'RelTol',acc1,'AbsTol',acc2);
    II2=integral(@(R) (K/S0).^(1-alpha-1i*R*exp(-1i*phi_star))./((-R*exp(-1i*phi_star)+1i*(alpha-1)).*(-R*exp(-1i*phi_star)+1i*alpha)).*FT_prob_dens(R*exp(-1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(-1i*phi_star),0,N/3,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha-1i*R*exp(-1i*phi_star))./((-R*exp(-1i*phi_star)+1i*(alpha-1)).*(-R*exp(-1i*phi_star)+1i*alpha)).*FT_prob_dens(R*exp(-1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(-1i*phi_star),N/3, N,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha-1i*R*exp(-1i*phi_star))./((-R*exp(-1i*phi_star)+1i*(alpha-1)).*(-R*exp(-1i*phi_star)+1i*alpha)).*FT_prob_dens(R*exp(-1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(-1i*phi_star),N,5*N,'RelTol',acc1,'AbsTol',acc2)+integral(@(R) (K/S0).^(1-alpha-1i*R*exp(-1i*phi_star))./((-R*exp(-1i*phi_star)+1i*(alpha-1)).*(-R*exp(-1i*phi_star)+1i*alpha)).*FT_prob_dens(R*exp(-1i*phi_star)-1i*alpha, T, S0, r, q, v0, kappa, theta, sigma, rho)*exp(-1i*phi_star),5*N,10*N,'RelTol',acc1,'AbsTol',acc2);
    CallPrice_v3(j)=-S0*exp(-r*T)/(2*pi)*real(II1+II2);
    dif(j)=CallPrice_v4(j)-CallPrice_v3(j);  
    
    IV1(j)=blsimpv(S0,K,r,T,CallPrice_v1(j)); %,'Method','jackel2016'
    IV2(j)=blsimpv(S0,K,r,T,CallPrice_v2(j)); %,'Method','jackel2016'
    IV3(j)=blsimpv(S0,K,r,T,CallPrice_v3(j),'Limit', 10); %,'Method','jackel2016'
    IV4(j)=blsimpv(S0,K,r,T,CallPrice_v4(j),'Limit', 10); %,'Method','jackel2016'
end 

% Top plot
nexttile
plot(KK,[IV1; IV2; IV3; IV4])
legend('IV1','IV2','IV3','IV4')
axis([0 K 0 5])
title('implied volatility')

% next plot
nexttile
plot(KK,IV1)
legend('IV1')
axis([0 K 0 5])
title('IV1')

% next plot
nexttile
plot(KK,IV2)
legend('IV2')
axis([0 K 0 5])
title('IV2')

% next plot
nexttile
plot(KK,IV3)
legend('IV3')
axis([0 K 0 5])
title('IV3')

% bottom plot
nexttile
plot(KK,IV4)
legend('IV4')
axis([0 K 0 5])
title('IV4')

