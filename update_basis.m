function B = update_basis(B,corr,C,W,D,Y,lamb,lambda,lambda_1,lambda_2,lambda_3)
%%run proximal gradient descent using current estimate

num_iter_max = 100; % max iter 
t = 0.0001; %lr for prox
err_inner = zeros(num_iter_max,1); %initialize

for iter = 1:num_iter_max
  
  %pre-initialize
  DG = zeros(size(B));
  
  %update gradient direction
  for j = 1:size(corr,1)
      
      %input variables
      Corr_j = reshape(corr(j,:,:),[size(corr,2),size(corr,3)]);
      D_k = reshape(D(j,:,:),[size(D,2),size(D,3)]);
      lamb_j =reshape(lamb(j,:,:),[size(lamb,2),size(lamb,3)]);
  
      %gradient update
      DG = DG -2*Corr_j*D_k +2*B*(D_k'*D_k) -D_k*diag(C(:,j)) +B*diag(C(:,j))*diag(C(:,j))- lamb_j*diag(C(:,j));
  
  end
  
  %proximal update    
  X_mat = B - t*DG/lambda_1;
  B = sign(X_mat).*(max(abs(X_mat)-t,0));

  %plot
  err_inner(iter) = error_compute(corr,B,C,Y,W,D,lamb,lambda,lambda_1,lambda_2,lambda_3);
  fprintf(' At B iteration %d || Error: %f \n',iter,err_inner(iter))   
  plot(1:iter,err_inner(1:iter),'b');
  hold on;
  drawnow;
  
  %check exit condition
  if(iter ==1)
   
      DG_init = DG;
  
  end
  
  if ((iter>5)&&((max(0,norm(DG,2)/norm(DG_init,2))< 10e-03)||(err_inner(iter)>err_inner(iter-5))))
      
      if(err_inner(iter)>err_inner(iter-5))
          fprintf((' Exiting due to increase in function value, adjust learning rate \n'))
      end
      
      break;
 
  end
  
end 