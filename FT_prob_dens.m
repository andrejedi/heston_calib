function ft=FT_prob_dens(w, t, S0,  mi, q, v0, kappa, theta, sigma, rho)
d=sqrt((rho*sigma*w*1i-kappa).^2+(w*1i+w.^2)*sigma^2);
g2=(kappa-rho*sigma*w*1i-d)./(kappa-rho*sigma*w*1i+d);
s1=exp(1i*w*(mi-q)*t);
s2=exp(theta*kappa/sigma^2*((kappa-rho*sigma*w*1i-d)*t-2*log((1-g2.*exp(-d*t))./(1-g2))));
s3=exp(v0/sigma^2*(kappa-rho*sigma*w*1i-d).*(1-exp(-d*t))./(1-g2.*exp(-d*t)));
ft=s1.*s2.*s3;
end
