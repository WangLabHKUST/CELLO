function c = mycolorgrad3(n,c1,c2,c3)
% Generate a color gradient map from c1 > c2 > c3

if mod(n,2) ~= 1
    error('Input n must be an odd number!')
end

m = median(1:n);

r = [linspace(c1(1),c2(1),m)';linspace(c2(1),c3(1),m)'];
g = [linspace(c1(2),c2(2),m)';linspace(c2(2),c3(2),m)'];
b = [linspace(c1(3),c2(3),m)';linspace(c2(3),c3(3),m)'];

c = [r,g,b];
c(m,:) = [];