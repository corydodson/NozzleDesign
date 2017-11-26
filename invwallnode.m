function [xOut,yOut,uOut,vOut] = invwallnode(t0,species,moleFracs,x,y,u,v,theta,d)

h0 = mixprop('h',species,moleFracs,t0);
r = mixprop('r',species,moleFracs)*1000;

u(3:4) = mean(u(1:2));
v(3:4) = mean(v(1:2));
dydx = diff(y(1:2))/diff(x(1:2));

N = 20;
err = 1e-6;
for i = 1:N
    u(3) = (u(3) + u(4))/2;
    v(3) = (v(3) + v(4))/2;
    theta(1:3) = atand(v(1:3)./u(1:3));
    vMag = sqrt(u.^2 + v.^2);
    
    h = h0 - vMag.^2/2000;
    t = tempfromprop(species,moleFracs,'h',h);
    g = mixprop('gamma',species,moleFracs,t);
    
    a = sqrt(r*g.*t);
    mu = asind(a./vMag);
    lp = tand(theta + mu);
    
    A = [-lp(3),1;-dydx,1];
    B = [y(4) - lp(3)*x(4);y(1) - dydx*x(1)];
    X = A\B;
    
    x(3) = X(1);
    y(3) = X(2);
    
    d12 = sqrt(diff(x(1:2))^2 + diff(y(1:2))^2);
    d13 = sqrt((x(3) - x(1))^2 + (y(3) - y(1))^2);
    
    u(3) = u(1) + diff(u(1:2))/d12*d13;
    v(3) = v(1) + diff(v(1:2))/d12*d13;
    
    theta(3) = atand(v(3)/u(3));
    vMag = sqrt(u.^2 + v.^2);

    h = h0 - vMag.^2/2000;
    t = tempfromprop(species,moleFracs,'h',h);
    g = mixprop('gamma',species,moleFracs,t);
    
    a = sqrt(r*g.*t);
    mu = asind(a./vMag);
    lp = tand(theta + mu);
    
    qp = u(3)^2 - a(3)^2;
    rp = 2*u(3)*v(3) - qp*lp(3);
    sp = d*a(3)^2*v(3)/y(3);
    tp = sp*diff(x(3:4)) + qp*u(3) + rp*v(3);
    
    A = [tand(theta(4)),-1;qp,rp];
    B = [0;tp];
    X = A\B;
    
    % check for convergence
    switch i ~= 1
        case 1
            notConv = abs([u(4)/X(1),v(4)/X(2)] - 1) > err;
            switch sum(sum(notConv)) == 0
                case 1
                    break
            end
    end
    
    u(4) = X(1);
    v(4) = X(2);
    
    % check for imaginary parts
    switch imag(u(4)) ~= 0 | imag(v(4)) ~= 0
        case 1
            return
    end
    
end

switch y(4) < y(2) + lp(2)*(x(4) - x(2))
    case 1
        return
end

xOut = x(4);
yOut = y(4);
uOut = X(1);
vOut = X(2);

end