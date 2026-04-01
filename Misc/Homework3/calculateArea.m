function [TriArea] = calculateArea(X1, X2, X3, Y1, Y2, Y3);

S1 = ([(X2-X1), (Y2-Y1), 0]);
S2 = ([(X3-X2), (Y3-Y2), 0]);

TriArea = norm(cross(S1, S2))/2;