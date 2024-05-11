function F = Model_ReturnFn(aprime,a,z,K_to_L,alpha,delta,crra,lam_hsv,tau_hsv)
% The return function is essentially the combination of the utility
% function and the constraints.

[r,w] = fun_prices(K_to_L,alpha,delta);

income = w*z+r*a;
taxes  = fun_hsv(income,lam_hsv,tau_hsv);
cons   = income+a-aprime-taxes; % Budget Constraint

F = -inf;

if cons>0
    % WARNING: this will not work if crra=1 and/or nu=1
    F = (cons^(1-crra)-1)/(1-crra);
end

end %end function