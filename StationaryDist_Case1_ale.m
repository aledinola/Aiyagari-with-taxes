function [mu_new] = StationaryDist_Case1_ale(Policy,n_d,n_a,n_z,pi_z, simoptions)

verbose = simoptions.verbose;
tol = 1e-8;

Policy1 = squeeze(gather(Policy)); % (n_a,n_z)

% Initial guess for mu
mu_old = ones(n_a,n_z)/(n_a*n_z); 

iter = 1;
err = tol+1;

while err>tol && iter<=10000

%mu_new = zeros(n_a,n_z);
% for z_c=1:n_z
%     for a_c=1:n_a
%         ap_opt = Policy1(a_c,z_c);
%         mu_new(ap_opt,:) = mu_new(ap_opt,:)+pi_z(z_c,:)*mu_old(a_c,z_c);
%     end
% end

mu_hat = zeros(n_a,n_z);
for z_c=1:n_z
    for a_c=1:n_a
        ap_opt = Policy1(a_c,z_c);
        mu_hat(ap_opt,z_c) = mu_hat(ap_opt,z_c)+mu_old(a_c,z_c);
    end
end
mu_new = mu_hat*pi_z;

% Compute distance
err = max(abs(mu_new-mu_old),[],'all');

if verbose==1
    fprintf('Iter = %d, err = %f \n',iter,err)
end

% Update
mu_old = mu_new;
iter = iter+1;


end %end while loop

end %end function