function [y, dy, ddy, f] = dmp_step( z, y, dy, ddy )
    global ay
    global by
    global c
    global h
    global n_dofs
    global goal
    global y_0
    global dt
    global w
    global tau
    global error_coupling
    

    % Generate Basis Function Activation %
    psi = exp(-h.*(z-c).^2);    %[n_bfs]
%     psi = exp(-h.*(mod(z,2*pi)-c).^2);   % 논문에 나온 식 대로 써봄.
%     psi = exp(h.*cos(z-c)-1); % 왜 일케하면 안되지??
    
    for n=1:n_dofs
        % generate forcing term
        front_term = z*(goal(n)-y_0(n));
        f(n) = front_term*(w(n, :)*psi')/sum(psi);
        % dmp acceleration
%         ddy(n) = (ay(n)*(by(n)*(goal(n)-y(n))-dy(n)/tau) + f(n))*tau; % ?
        ddy(n) = (ay(n)*(by(n)*(goal(n)-y(n))-dy(n)) + f(n))*(tau^2);
%         dy(n) = dy(n) + ddy(n)*tau*dt*error_coupling;
        dy(n) = dy(n) + ddy(n)*dt;
        y(n) = y(n) + dy(n)*dt;
    end
    
end