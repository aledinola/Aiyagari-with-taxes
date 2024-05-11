function taxes = fun_hsv(income,lambda,tau)
% Income taxes, modelled following Heathcote, Storesletten and Violante

income_pos = max(income,0);
taxes  = income_pos-lambda*income_pos^(1-tau);

end