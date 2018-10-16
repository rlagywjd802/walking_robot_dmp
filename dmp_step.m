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
    
    for n=1:n_dofs
        % generate forcing term
        front_term = z*(goal(n)-y_0(n));
        f(n) = front_term*(w(n, :)*psi')/sum(psi);
        % dmp acceleration
        ddy(n) = (ay(n)*(by(n)*(goal(n)-y(n))-dy(n)/tau) + f(n))*tau;
        dy(n) = dy(n) + ddy(n)*tau*dt*error_coupling;
        y(n) = y(n) + dy(n)*dt*error_coupling;
    end
    
end

