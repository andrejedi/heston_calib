function [ sigmas ] = impliedvol_c(MarketPrice,S0,K,T,r)

    MarketPrice = MarketPrice(:);
    S0 = S0(:);
    K = K(:);
    T = T(:);
    sigmas = zeros(size(K,1),1);

    for i=1:size(K,1)
        implvolbs = blsimpv(S0(i),K(i),r,T(i),MarketPrice(i));
        sigmas(i) = implvolbs;
    end
end