function [q1,q2]=IK(x2, y2)
    global l1;
    global l2;
    
    c2 = (x2^2+y2^2-l1^2-l2^2)/(2*l1*l2) ;
    s2 = sqrt(1-c2^2);
    q2 = atan2(s2,c2);
   
%    D = det([l1+l2*c2, -l2*s2; l2*s2, l1+l2*c2]);
%    D1 = det([x2, -l2*s2; y2, l1+l2*c2]);
%    D2 = det([l1+l2*c2, x2; l2*s2, y2]);

%    c1 = D1/D;
%    s1 = D2/D;
    
%    q1 = atan2(s1,c1);

    q1 = atan2(y2,x2) - atan2(l2*s2, l1+l2*c2);

end
