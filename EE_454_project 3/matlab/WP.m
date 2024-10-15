function [Win] = WP(p)
Win = [
1+p(1)  p(3)    p(5);
p(2)    1+p(4)  p(6);
0       0       1   ];
end