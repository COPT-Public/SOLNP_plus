function solnp_tesths()
clear 
persistent str
op.tol = 1e-3;
op.tol_con = 1e-3; 
op.ls_time = 5;
op.min_iter = 3;
op.c_i = 5;op.k_i = 3;
op.c_r = 1;op.k_r = 3;


% op.re_time = 0;op.ls_time = 1000;op.noise = 0; % The origin SOLNP setting
rep = 50;
data = [];

str = pwd;
str = strcat(str,'\');  
totalt = 0;

s = strcat(str,'P11\');
cd(s);
disp("P11")
[f,constraint, count_cost,t,s] = P11(op,rep);
totalt = t + totalt;
data = [data ; f,constraint, count_cost,t,s];

s = strcat(str,'P26\');
cd(s);
disp("P26")
[f,constraint, count_cost,t,s] =P26(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P28\');
cd(s);
disp("P28")
[f,constraint, count_cost,t,s] =P28(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P38\');
cd(s);
disp("P38")
[f,constraint, count_cost,t,s] =P38(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P40\');
cd(s);
disp("P40")
[f,constraint, count_cost,t,s] =P40(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P46\');
cd(s);
disp("P46")
[f,constraint, count_cost,t,s] =P46(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P56\');
cd(s);
disp("P56")
[f,constraint, count_cost,t,s] =P56(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P78\');
cd(s);
disp("P78")
[f,constraint, count_cost,t,s] =P78(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P79\');
cd(s);
disp("P79")
[f,constraint, count_cost,t,s] =P79(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P80\');
cd(s);
disp("P80")
[f,constraint, count_cost,t,s] =P80(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P81\');
cd(s);
disp("P81")
[f,constraint, count_cost,t,s] =P81(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P83\');
cd(s);
disp("P83")
[f,constraint, count_cost,t,s] =P83(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P84\');
cd(s);
disp("P84")
[f,constraint, count_cost,t,s] =P84(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P93\');
cd(s);
disp("P93")
[f,constraint, count_cost,t,s] =P93(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

s = strcat(str,'P106\');
cd(s);
disp("P106")
[f,constraint, count_cost,t,s] =P106(op,rep);
totalt = t + totalt;
data = [data;f,constraint, count_cost,t,s];

cd(str);
save("noise-data.mat","data",'-mat');
fprintf("total time = %e",totalt);
end