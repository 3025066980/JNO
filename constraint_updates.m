function [D_upd,lamb_upd]= constraint_updates(corr,B_upd,C_upd,lamb,lr)
%update constraint variables D and lamb

num_iter_max=100; %max iter

parfor k= 1:size(C_upd,2) %decoupled patient wise

     %input variables
     Corr_k = reshape(corr(k,:,:),[size(corr,2),size(corr,3)]);
     lamb_k =reshape(lamb(k,:,:),[size(lamb,2),size(lamb,3)]);
     
     for c=1:num_iter_max
        
        %primal dual updates
        D_k = (B_upd*diag(C_upd(:,k))+ 2*Corr_k*B_upd - lamb_k)*pinv(eye(size(B_upd'*B_upd))+2*(B_upd'*B_upd));
        lamb_k = lamb_k + (0.5^(c-1))*lr*(D_k - B_upd*diag(C_upd(:,k)));
        
        %check exit condition
        if (c ==1)
            grad_norm_init = norm(D_k - B_upd*diag(C_upd(:,k)),2);
        end
        
        if (norm(D_k - B_upd*diag(C_upd(:,k)),2)/grad_norm_init<10e-03)
           break;
        end
        
     end
     
     %update D and lamb
     lamb_upd(k,:,:)= lamb_k;
     D_upd(k,:,:) =D_k;
     
end

end