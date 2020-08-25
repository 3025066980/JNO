function [B,C,W,D,lamb] = altmin_run(corr,B_init,C_init,W_init,D_init,Y,lamb_init,lambda,lambda_1,lambda_2,lambda_3,lr)
%%runs alternating minimization with initial solution

%Initilise
num_iter = 100;
B_old = B_init;
C_old = C_init;
D_old = D_init;
W_old = W_init;
lamb_old = lamb_init;

%thresh for exit
thresh = 10e-04;

%plot   
err_out = zeros(num_iter,1);
plot(0,0)
title('Gradent Descent Run ')
xlabel('Number of iterations')
ylabel('Value of obejctive function')

%Iterate
    for i = 1:num_iter
        
        err_out(i)= error_compute(corr,B_old,C_old,Y,W_old,D_old,lamb_old,lambda,lambda_1,lambda_2,lambda_3);
      
        fprintf(' At iteration %d || Error: %f \n',i,err_out(i))
        plot(1:i,err_out(1:i),'r');
        hold on;
        drawnow;
            
        [B,C,D,W,lamb] = altmin_step(corr,B_old,C_old,W_old,D_old,lamb_old,Y,lambda,lambda_1,lambda_2,lambda_3,lr); 
        
        %update
        B_old = B;
        C_old = C;
        D_old = D;
        W_old = W;
        lamb_old = lamb;
    
        %exit condition
        if (i>1 && (abs((err_out(i)-err_out(i-1))) < thresh || (err_out(i)-err_out(i-1)>10)))
           
            if (err_out(i)>err_out(i-1))
               fprintf('Exiting due to increase in function value, adjust learning rate')
            end
            
            break;
        
        end
    
    end
    
    
end