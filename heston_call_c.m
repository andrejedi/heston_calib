function [ CallPrice_v1 ] = heston_call_c( S0,K,T,q,alpha,r,v0,kappa,theta,sigma,rho,N )
    % returns Heston Price for vanilla european options
    %
    % S: underlying spot price
    % K: strike price
    % T: date at end of maturity
    % t: start date (T-t is time to maturity)
    % r: interest rate
    % q: dividend rate
    % sigma: volatility of returns of underlying
    
    S0=S0(:);
    K=K(:);
    T=T(:);

    CallPrice_v1=zeros(size(K,1),1);
    
    for i = 1:numel(K)
        I1=integral(@(w) -(K(i)/S0(i)).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T(i), S0(i), r, q, v0, kappa, theta, sigma, rho),0,N)+...
        integral(@(w) -(K(i)/S0(i)).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T(i), S0(i), r, q, v0, kappa, theta, sigma, rho),N,5*N)+...
        integral(@(w) -(K(i)/S0(i)).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T(i), S0(i), r, q, v0, kappa, theta, sigma, rho),5*N,10*N);    
        I2=integral(@(w) -(K(i)/S0(i)).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T(i), S0(i), r, q, v0, kappa, theta, sigma, rho),-N,0)+...
        integral(@(w) -(K(i)/S0(i)).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T(i), S0(i), r, q, v0, kappa, theta, sigma, rho),-5*N,-N)+...
        integral(@(w) -(K(i)/S0(i)).^(1-alpha+1i*w)./((w+1i*(alpha-1)).*(w+1i*alpha)).*FT_prob_dens(-w-1i*alpha, T(i), S0(i), r, q, v0, kappa, theta, sigma, rho),-10*N,-5*N);
        CallPrice_v1(i)=S0(i)*exp(-r*T(i))/(2*pi)*real(I1+I2); 
    end
end