function [theta]=findAngle(a,b)

% input is two 3x1 vectors
% output is the angle between those vectors

if numel(a)==3
    aDotb=a(1)*b(1)+a(2)*b(2)+a(3)*b(3);
    anorm=sqrt(a(1)^2+a(2)^2+a(3)^2);
    bnorm=sqrt(b(1)^2+b(2)^2+b(3)^2);
else
    aDotb=a(1)*b(1)+a(2)*b(2);
    anorm=sqrt(a(1)^2+a(2)^2);
    bnorm=sqrt(b(1)^2+b(2)^2);
end
    

theta=acos(aDotb/(anorm*bnorm));

