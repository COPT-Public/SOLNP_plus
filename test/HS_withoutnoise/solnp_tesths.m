function solnp_tesths()
clear 
persistent str
op.tol = 1e-4;
op.tol_con = 1e-4; 
op.ls_time = 5;
op.min_iter = 10;
% op.c_i = 30;op.k_i = 3;
% op.c_r = 10;op.k_r = 9;

rep = 50;
data = [];

str = pwd;
str = strcat(str,'\');  
totalt = 0;

s = strcat(str,'P11\');
cd(s);
disp("P11")
[f,constraint, count_cost,t] = P11(op,rep);
totalt = t + totalt;
data = [data ; f,constraint, count_cost,t];

s = strcat(str,'P26\');
cd(s);
disp("P26")
[f,constraint, count_cost,t] =P26(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P28\');
cd(s);
disp("P28")
[f,constraint, count_cost,t] =P28(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P38\');
cd(s);
disp("P38")
[f,constraint, count_cost,t] =P38(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P40\');
cd(s);
disp("P40")
[f,constraint, count_cost,t] =P40(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P46\');
cd(s);
disp("P46")
[f,constraint, count_cost,t] =P46(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P56\');
cd(s);
disp("P56")
[f,constraint, count_cost,t] =P56(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P78\');
cd(s);
disp("P78")
[f,constraint, count_cost,t] =P78(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P79\');
cd(s);
disp("P79")
[f,constraint, count_cost,t] =P79(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P80\');
cd(s);
disp("P80")
[f,constraint, count_cost,t] =P80(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P81\');
cd(s);
disp("P81")
[f,constraint, count_cost,t] =P81(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P83\');
cd(s);
disp("P83")
[f,constraint, count_cost,t] =P83(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P84\');
cd(s);
disp("P84")
[f,constraint, count_cost,t] =P84(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P93\');
cd(s);
disp("P93")
[f,constraint, count_cost,t] =P93(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

s = strcat(str,'P106\');
cd(s);
disp("P106")
[f,constraint, count_cost,t] =P106(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t];

cd(str)
save("withoutnoise-data.mat","data",'-mat');
fprintf("total time = %e",totalt);
end