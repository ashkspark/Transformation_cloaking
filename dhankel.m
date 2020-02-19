function [h] = dhankel(nu,z)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

h=nu*besselh(nu,z)./z-besselh(nu+1,z);

end

