function [ f1 ] = Centhalpy(  h , x ,T, P  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[~,~,h1,~,~]= srk(x,T,P,1);
f1 = h - h1;
end
