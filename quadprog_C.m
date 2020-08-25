function C_upd = quadprog_C(B_upd,corr,C,W,Y,D,lamb,lambda,lambda_2)
%run quadratic program to update the coefficients

%pre initialize
C_upd = zeros(size(C));

%decoupling across patients
parfor m = 1:size(corr,1)
    
    %quad term
    H = diag(diag(B_upd'*B_upd)) + 2*(lambda*(W*W')+ 2*lambda_2* eye(size(B_upd'*B_upd)));
    
    D_m = reshape(D(m,:,:),[size(D,2),size(D,3)]);
    lamb_m =reshape(lamb(m,:,:),[size(lamb,2),size(lamb,3)]);
    
    L1 = D_m'*B_upd;
    L2 = lamb_m'*B_upd;
    
    %linear term
    f = -diag(L1)-diag(L2)-2*lambda*Y(m)*W;
    
    %constraints- non neg
    A = -eye(size(C,1));
    b = zeros(size(C,1),1);
    
    c_m = quadprog(H,f,A,b);
    C_upd(:,m) = c_m;
    
end

end