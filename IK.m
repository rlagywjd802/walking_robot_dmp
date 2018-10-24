function [q1,q2]=IK(x2, y2)
    global l1;
    global l2;
    
    % Siciliano 92pg
    c2 = (x2^2+y2^2-l1^2-l2^2)/(2*l1*l2);
    s2 = -sqrt(1-c2^2); % elbow-up
    q2 = atan2(s2,c2);
    
    c1 = ((l1+l2*c2)*y2+l2*s2*x2)/(x2^2+y2^2);
    s1 = ((l1+l2*c2)*y2-l2*s2*x2)/(x2^2+y2^2);
    q1 = atan2(s1, c1);

    q2 = -q2;
    q1 = q1 - (3./2.)*pi;
    
    
    

end
