function psi0 = GetPsi0(x2T, k, L, alpha, T)
x1_T = L + alpha*T/(2*(1+k)) + alpha*exp(-(1+k)*T)/(2*(1+k)^2) - ...
        alpha/(2*(1+k)^2);
    
%x2T = (S - eps) or x2T = (S - eps)
x2_T = x2T + (alpha/(2*(1+k))) * (1-exp(-(1+k)*T));
    
    c1 = T/(1+k)^2 - exp((1+k)*T)/(2*(1+k)^3) + ...
        exp(-(1+k)*T)/(2*(1+k)^3); %the coefficient at psi_1^0 at x1(T)
    d1 = exp((1+k)*T)/(2*(1+k)^2) + exp(-(1+k)*T)/(2*(1+k)^2) - ...
        1/((1+k)^2); %the coefficient at psi_2^0 at x1(T)
    
    c2 = 1/(1+k)^2 - exp((1+k)*T)/(2*(1+k)^2) - ...
        exp(-(1+k)*T)/(2*(1+k)^2); %the coefficient at psi_1^0 at x2(T)
    d2 = exp((1+k)*T)/(2*(1+k)) - ...
        exp(-(1+k)*T)/(2*(1+k)); %the coefficient at psi_2^0 at x2(T)
    
    A1 = [c1 d1; c2 d2];
    B1 = [x1_T; x2_T]; 
    
    if det(A1) ~= 0
        %psi0 = inv(A1)*B1; % vector psi(0)
        psi0 = A1\B1; % vector psi(0)
    else
        error('det(a) = 0: var 1, psi_2(T) > 0')
    end
end