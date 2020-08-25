function [B_upd,C_upd,D_upd,W_upd,lamb_upd] = altmin_step(corr,B,C,W,D,lamb,Y,lambda,lambda_1,lambda_2,lambda_3,lr)
%%Given the current values of the iterates, performs a single step of alternating minimisation

%% B update
fprintf('Optimise B \n')

%prox updates
B_upd = update_basis(B,corr,C,W,D,Y,lamb,lambda,lambda_1,lambda_2,lambda_3);

fprintf(' At final B iteration || Error: %f \n',error_compute(corr,B_upd,C,Y,W,D,lamb,lambda,lambda_1,lambda_2,lambda_3))   

%% C update

fprintf('Optimise C \n')

%quad prog solution
C_upd = quadprog_C(B_upd,corr,C,W,Y,D,lamb,lambda,lambda_2);

fprintf(' Step C || Error: %f \n',error_compute(corr,B_upd,C_upd,Y,W,D,lamb,lambda,lambda_1,lambda_2,lambda_3));


%% W update

fprintf('Optimise W \n')

%closed form solution for regression weights

if  (lambda~=0) %regular
W_upd = ((C_upd*C_upd')+(lambda_3/lambda)*eye(size(C*C')))\(C_upd*Y);
else % without regression term
W_upd = zeros(size(B_upd,2),1);
end

fprintf(' Step W || Error: %f \n',error_compute(corr,B_upd,C_upd,Y,W_upd,D,lamb,lambda,lambda_1,lambda_2,lambda_3));

%% constraint updates

fprintf('Optimise D and lambda \n')

%constraint updates
[D_upd,lamb_upd] = constraint_updates(corr,B_upd,C_upd,lamb,lr);

fprintf(' Step D/lamb || Error: %f \n',error_compute(corr,B_upd,C_upd,Y,W_upd,D_upd,lamb_upd,lambda,lambda_1,lambda_2,lambda_3));
       
end
     

   
    