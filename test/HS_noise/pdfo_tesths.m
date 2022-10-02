function pdfo_tesths()
clear 
persistent str


str = pwd;
str = strcat(str,'\');  
totalt = 0;

s = strcat(str,'P11\');
cd(s);
disp("PDFOP11")
eval(" PDFOP11;");
%totalt = t + totalt;


s = strcat(str,'P26\');
cd(s);
disp("PDFOP26")
eval("PDFOP26;");
%totalt = t + totalt;


% s = strcat(str,'P28\');
% cd(s);
% disp("PDFOP28")
% eval("PDFOP28;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P38\');
% cd(s);
% disp("PDFOP38")
% eval("PDFOP38;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P40\');
% cd(s);
% disp("PDFOP40")
% eval("PDFOP40;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P46\');
% cd(s);
% disp("PDFOP46")
% eval("PDFOP46;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P56\');
% cd(s);
% disp("PDFOP56")
% eval("PDFOP56;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P78\');
% cd(s);
% disp("PDFOP78")
% eval("PDFOP78;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P79\');
% cd(s);
% disp("PDFOP79")
% eval("PDFOP79;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P80\');
% cd(s);
% disp("PDFOP80")
% eval("PDFOP80;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P81\');
% cd(s);
% disp("PDFOP81")
% eval("PDFOP81;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P83\');
% cd(s);
% disp("PDFOP83")
% eval("PDFOP83;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P84\');
% cd(s);
% disp("PDFOP84")
% eval("PDFOP84;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P93\');
% cd(s);
% disp("PDFOP93")
% eval("PDFOP93;");
% %totalt = t + totalt;
% 
% 
% s = strcat(str,'P106\');
% cd(s);
% disp("PDFOP106")
% eval("PDFOP106;");
% %totalt = t + totalt;
% 
% 
% cd(str);
% % save("noise-data-pdfo.mat","data",'-mat');
fprintf("total time = %e",totalt);
end