function [y]  =LV(t,x,r,K,A,Delta) 
%implementation of Eqs. (2) and (3)
y=-Delta.*x+x.*(r-r.*x/K) + x.*(A*x);
end
