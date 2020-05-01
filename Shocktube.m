clear
clc
close all

input   = loadinVar();

s1      = EulerFiniteVolume1D(input);

s1.solveExplicit();
