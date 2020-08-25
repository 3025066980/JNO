function err = error_compute(corr,B,C,Y,W,D,lamb,lambda,lambda_1,lambda_2,lambda_3)
%%Computes the error at the current interation given the values of the iterates at the instant

%initialise
fit_err = 0;
const_err =0 ;
aug_lag_err =0;

for n = 1:size(corr,1)
    
    D_n = reshape(D(n,:,:),[size(D,2),size(D,3)]);
    lamb_n = reshape(lamb(n,:,:),[size(lamb,2),size(lamb,3)]);
    
    X = reshape((corr(n,:,:)),[size(corr,2),size(corr,3)]) -D_n*B'; 
    fit_err = fit_err + norm(X,'fro').^2;
    
    const_err = const_err + trace(lamb_n'*(D_n-B*diag(C(:,n))));
    aug_lag_err = aug_lag_err +0.5*norm((D_n-B*diag(C(:,n))),'fro').^2;
 
end

err = fit_err + const_err + aug_lag_err + lambda_1* norm(B,1) + lambda* norm(C'*W-Y,2).^2 + lambda_2* norm(C,'fro').^2+ lambda_3*norm(W,'fro').^2;

end