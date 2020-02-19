function [l] = dbessely(nu,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
l=0.5.*(bessely(nu-1,z) - bessely(nu+1,z));

end

