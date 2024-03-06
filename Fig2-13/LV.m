function [y]  =LV(t,x,r,K,A,Delta) 
y=-Delta.*x+x.*(r-r.*x/K) + x.*(A*x);
%y=-Delta.*x+x.*r+ x.*(A*x);
end

