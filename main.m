clear all; close all

%% parameters
lr = 0.0001; 
lambda = 1; %regression
lambda_1 = 30;    %sparsity
lambda_2 = 0.2;   % regularizer C
lambda_3 = 1;     % ridge regularizer W
net = 8;     %networks

m= 116; %AAL atlas ROIs

%% Initialising parameters

tic;

%load inputs
foldername = '/home/niharika-shimona/Documents/Projects/Autism_Network/JNO/Data/ADOS/';
filename = strcat(foldername,'/data.mat');   
load(filename)
    
%initialize
corr_mean = reshape(mean(corr_train,1),[m,m]); %average

[V,D] = eig(corr_mean); 
[~,permutation] = sort(diag(D),'descend');
D=D(permutation,permutation);V=V(:,permutation);

B_init = V(:,1:net);%init basis
C_init = quadestimate_Ctest(B_init,corr_train,lambda_2); %init coeffs
W_init = randn(net,size(Y_train,2)); %init regression weights

%init constraint variables
[D_init,lamb_init] = init_const(B_init,C_init) ;

%print params
fprintf('\n sparsity penalty : %f; networks: %d \n',lambda_1,net)
fprintf('\n tradeoff penalty: %f, C penalty : %f; W penalty: %f \n',lambda,lambda_2,lambda_3)
 
%run alt min procedure
[B_gd,C_gd,W_gd,D_gd,lamb_gd] = altmin_run(corr_train,B_init,C_init,W_init,D_init,Y_train,lamb_init,lambda,lambda_1,lambda_2,lambda_3,lr);

%test coeffs
C_gd_train = quadestimate_Ctest(B_gd,corr_train,lambda_2);
C_gd_test = quadestimate_Ctest(B_gd,corr_test,lambda_2);

%testing performance
Y_train_meas = Y_train;
Y_train_pred = C_gd'*W_gd;
Y_train_pred_reg = C_gd_train'*W_gd;

Y_test_meas = Y_test;
Y_test_pred = C_gd_test'*W_gd;

%save outputs
output_filename =  strcat(foldername,'/Outputs/models.mat');
save(output_filename,'B_gd','C_gd','C_gd_test','C_gd_train','D_gd','lamb_gd');

output_perf =  strcat(foldername,'/Outputs/Performance.mat');
save(output_perf,'Y_train_meas','Y_train_pred','Y_train_pred_reg','Y_test_meas','Y_test_pred');
