function C_upd = quadestimate_Ctest(B,corr,lambda_2)
%calculate coefficients at testing time

%pre-initialize
C_upd = zeros(size(B,2),size(corr,1));

parfor m = 1:size(corr,1)
   
    %quad term
    H = 2*((B'*B).^2 + lambda_2* eye(size(B'*B)));
   
    %linear term
    Corr_mat = reshape(corr(m,:,:),[size(corr,2),size(corr,3)]);
    M = -2*(B'*(Corr_mat)*B);
    f = diag(M) ;
    
    %constraints- non neg
    A = -eye(size(B,2));
    b = zeros(size(B,2),1);
    c_m = quadprog(H,f,A,b);
    
    C_upd(:,m) = c_m;
    
end

end