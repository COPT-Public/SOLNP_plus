 function test()
persistent str
str = pwd;
str = strcat(str,'\');

s = strcat(str,'P11\');
cd(s);
disp("NOMADP11")
eval("NOMADP11");

s = strcat(str,'P26\');
cd(s);
disp("NOMADP26")
eval("NOMADP26");

s = strcat(str,'P38\');
cd(s);
disp("NOMADP38")
eval("NOMADP38");

s = strcat(str,'P40\');
cd(s);
disp("NOMADP40")
eval("NOMADP40");

s = strcat(str,'P46\');
cd(s);
disp("NOMADP46")
eval("NOMADP46");

s = strcat(str,'P56\');
cd(s);
disp("NOMADP56")
eval("NOMADP56");

s = strcat(str,'P78\');
cd(s);
disp("NOMADP78")
eval("NOMADP78");

s = strcat(str,'P79\');
cd(s);
disp("NOMADP79")
eval("NOMADP79");

s = strcat(str,'P80\');
cd(s);
disp("NOMADP80")
eval("NOMADP80");

s = strcat(str,'P81\');
cd(s);
disp("NOMADP81")
eval("NOMADP81");

s = strcat(str,'P83\');
cd(s);
disp("NOMADP83")
eval("NOMADP83");

s = strcat(str,'P84\');
cd(s);
disp("NOMADP84")
eval("NOMADP84");

s = strcat(str,'P93\');
cd(s);
disp("NOMADP93")
eval("NOMADP93");

s = strcat(str,'P106\');
cd(s);
disp("NOMADP106")
eval("NOMADP106");
end