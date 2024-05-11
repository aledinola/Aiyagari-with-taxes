function [r,w] = fun_prices(K_to_L,alpha,delta)

r = alpha*K_to_L^(alpha-1)-delta;
w = (1-alpha)*K_to_L^alpha;

end