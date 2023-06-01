clc;
clear all;
Func_num = input('Input the function number: ');
global initial_flag
initial_flag = 0;
Par     = Cal_par(Func_num);
D       = Par.n
g       = Par.g
h       = Par.h
xmin    = Par.xmin
xmax    = Par.xmax
Solution = (xmax-xmin).*rand(1,D)+xmin
[f,g,h] = cec20_func(Solution,Func_num)