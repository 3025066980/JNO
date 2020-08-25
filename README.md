   

# JNO

Code for JNO model introduced in https://www.sciencedirect.com/science/article/pii/S105381191930905X and https://link.springer.com/chapter/10.1007%2F978-3-030-00931-1_19

![JNO](https://github.com/Niharika-SD/JNO/blob/master/Images/Connectomics_and_Clinical_Severity.PNG)
##Instructions

1. Open MATLAB (needs parallel computing toolbox)
2. Run main.m


#DATA ORGANIZATION

Main directory contents:

1. main.m #main script
2. altmin_run.m #alternating minimization loop
3. altmin_step.m #runs a single step of alternating minimization
4. update_basis.m #proximal gradient descent update
5. quadprog_C.m #quadratic programming update
6. constraint_updates.m #augmented lagrangian updates
7. error_compute.m #computes JNO error for a single iterate
8. quadestimate_Ctest.m #quadratic programming at test
9. init_constraints.m #initialize constraint variables D and lamb

Folders:

~\Data\ADOS\data.mat
~\Data\ADOS\Outputs\Performance.mat
~\Data\ADOS\Outputs\models.mat

~\Data\SRS\data.mat
~\Data\SRS\Outputs\Performance.mat
~\Data\SRS\Outputs\models.mat

~\Data\Praxis\data.mat
~\Data\Praxis\Outputs\Performance.mat
~\Data\Praxis\Outputs\models.mat
   

