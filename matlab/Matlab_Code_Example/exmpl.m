%%%%%% RENAME CAFunction.exe_ TO CAFunction.exe!!!!!!!!

clear all;

cl_num = 10;
Eita = 2.5;
MaxIter = 50;

MaxPts =  zeros(50,1); %soft merging condition
MaxPts2 = [ones(25,1); 5*ones(25,1)]; %crisp merging condition
 
TAU = 10;
EXPINC = 25;
EXPDEC = 35;
m_p = 2;

f = [rand(20,2); rand(20,2)+1; rand(20,2)+2;]; % features  (example with 3 clusters), total of 60 feature vectors.
c = rand(cl_num, 2); %random initial centers 

[cl_num cl_centers cl_assigment progress_log] = CAalg(f, c, Eita, MaxIter, MaxPts, MaxPts2, TAU, EXPINC, EXPDEC, m_p);

cl_num
cl_centers
% progress_log
