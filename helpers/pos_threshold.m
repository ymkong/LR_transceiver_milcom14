function [ A ] = pos_threshold( A )
%POS_THRESHOLD Summary of this function goes here
%   Detailed explanation goes here

A(A < 0) = 0;
end

