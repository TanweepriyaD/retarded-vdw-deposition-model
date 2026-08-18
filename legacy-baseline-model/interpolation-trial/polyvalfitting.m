clc
clear all
close all

vdata =  xlsread('dfPsWtr_1M.xlsx','vdwdata_1M');
L =vdata(:,1)*1e-09;
vdwr =vdata(:,2);

coefficients = polyfit(L, vdwr, 2);

interpolatedValues = polyval(coefficients, xq); % xq vector of points where you want to interpolate
