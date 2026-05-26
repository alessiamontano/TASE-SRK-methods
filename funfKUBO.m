function yp = funfKUBO(y)
global a b
yp =  [-a*y(2)-b^2/2*y(1);a*y(1)-b^2/2*y(2)];
end