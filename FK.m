function [x] = FK(q)
    % Forward Kinematics of two link robot arm
    % q : joint angle in degree [N, 2]
    % x : EE of two link        [N, 2]    
    global R2D
    global D2R
    global l1
    global l2
    
    q1 = 3./2.*pi + q(:, 1)*D2R;
    q2 = - q(:, 2)*D2R;
    x1 = l1*cos(q1);
    y1 = l1*sin(q1);
    x2 = x1 + l2*cos(q1+q2);
    y2 = y1 + l2*sin(q1+q2);
    
%     q = q*D2R;
%     x1 = l1*sin(q(:,1));
%     y1 = -l1*cos(q(:,1));
%     x2 = x1 + l2*sin(q(:,1)-q(:,2));
%     y2 = y1 - l2*cos(q(:,1)-q(:,2));
    
    x = [x2, y2];
end

