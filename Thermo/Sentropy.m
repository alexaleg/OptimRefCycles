function [ f1 ] = Sentropy( S, x ,T, P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,~,~,~,S1]= srk(x,T,P,1);
f1 = S - S1;
end
