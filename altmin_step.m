function [B_upd,C_upd,K_upd,D_upd,lamb_upd] = altmin_step(corr,B,C,K,D,lamb,Y,lambda,lambda_1,lambda_2,lambda_3,lr,sigma,w_p,p)
%Given the current values of the iterates, performs a single step of alternating minimisation

%% B update
fprintf('Optimise B \n')

%prox updates
B_upd = update_basis(B,corr,C,W,D,Y,lamb,lambda,lambda_1,lambda_2,lambda_3);

fprintf(' At final B iteration || Error: %f \n',error_compute(corr,B_upd,C,C,K,Y,D,lamb,lambda,lambda_1,lambda_2,lambda_3,sigma,w_p,p,0))   

%% C update

fprintf('Optimise C \n')

%quad prog solution
alpha = pinv(K+(lambda_3/lambda)*eye(size(C,2)))*Y; %pre-compute alpha term
C_upd = coefficient_updates(corr,B_upd,C,Y,D,lamb,alpha,lambda,lambda_2,sigma,w_p,p);

fprintf(' Step C || Error: %f \n',error_compute(corr,B_upd,C_upd,C,K,Y,D,lamb,lambda,lambda_1,lambda_2,lambda_3,sigma,w_p,p,1));


%% Kernel Update- regression term

K_upd = compute_kernel(C_upd,sigma,w_p,p);

%check rank of matrix
fprintf(' Rank of Kernel Matrix: %d \n', rank(K_upd))

fprintf(' Step W-kernel || Error: %f \n',error_compute(corr,B_upd,C_upd,C_upd,K_upd,Y,D,lamb,lambda,lambda_1,lambda_2,lambda_3,sigma,w_p,p,0));

%% constraint updates

fprintf('Optimise D and lambda \n')

%constraint updates
[D_upd,lamb_upd] = constraint_updates(corr,B_upd,C_upd,lamb,lr);

fprintf(' Step D/lamb || Error: %f \n',error_compute(corr,B_upd,C_upd,C_upd,K_upd,Y,D_upd,lamb_upd,lambda,lambda_1,lambda_2,lambda_3,sigma,w_p,p,0));
       
end
     

   
    