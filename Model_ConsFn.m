function cons = Model_ConsFn(aprime,a,z,K_to_L,alpha,delta,lam_hsv,tau_hsv)
% The return function is essentially the combination of the utility
% function and the constraints.

[r,w] = fun_prices(K_to_L,alpha,delta);

income = w*z+r*a;
taxes  = fun_hsv(income,lam_hsv,tau_hsv);
cons   = income+a-aprime-taxes; % Budget Constraint

end %end function