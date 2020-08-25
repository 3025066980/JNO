function [D_init,lamb_init] = init_const(B_init,C_init)
%initialize the constraint variables

%pre-allocate
D_init = zeros(size(C_init,2),size(B_init,1),size(B_init,2));

%solution in feasible set
parfor j = 1:size(C_init,2)
    D_init(j,:,:) = B_init*diag(C_init(:,j));
end 

lamb_init = zeros(size(D_init));

end