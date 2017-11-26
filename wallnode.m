function [xOut,yOut,uOut,vOut] = wallnode(t0,species,moleFracs,x,y,u,v,d,f)

h0 = mixprop('h',species,moleFracs,t0);
r = mixprop('r',species,moleFracs)*1000;

L = length(x);

minus = logical([ones(L - 1,1);0]);
plus = logical([0;ones(L - 1,1)]);

xm = x(minus);
xmCalc = xm;
xp = x(plus);

ym = y(minus);
ymCalc = ym;
yp = y(plus);
ypCalc = yp;

um = u(minus);
umCalc = um;
up = u(plus);
upCalc = up;

vm = v(minus);
vmCalc = vm;
vp = v(plus);
vpCalc = vp;

N = 20;
err = 1e-4;
notConv = true(1,4);
for i = 1:N
    
    ypCalc = (yp + ypCalc)/2;
    umCalc = (um + umCalc)/2;
    upCalc = (up + upCalc)/2;
    vmCalc = (vm + vmCalc)/2;
    vpCalc = (vp + vpCalc)/2;
    
    uCalc = [umCalc;upCalc(end)];
    vCalc = [vmCalc;vpCalc(end)];
    
    vMag = sqrt(uCalc.^2 + vCalc.^2);
    h = h0 - vMag.^2/2000;
    t = tempfromprop(species,moleFracs,'h',h);
    g = mixprop('gamma',species,moleFracs,t);
    
    a = sqrt(r*g.*t);
    am = a(minus);
    
    mu = asind(a./vMag);
    mum = mu(minus);
    
    theta = atand(vCalc./uCalc);
    thetam = theta(minus);
    
    lm = tand(thetam - mum);
    
    q = uCalc.^2 - a.^2;
    qm = q(minus);
    
    rm = 2*umCalc.*vmCalc - qm.*lm;
    
    sm = d*am.^2.*vmCalc./ymCalc;
    
    % find x and y
    xmCalc = fzero(@(p) (polyval(f,p) - ym - lm*(p - xm)), xm);
    xpCalc = xmCalc;
    ymCalc = polyval(f,xmCalc);
    ypCalc = ymCalc;
    
    % find u and v
    A = [qm,rm;...
        tand(thetam),-1];
    B = [sm*(xmCalc - xm) + qm*um + rm*vm;...
        0];
    
    X = A\B;
    umCalc = X(1:2:end);
    upCalc = umCalc;
    vmCalc = X(2:2:end);
    vpCalc = vmCalc;
    
    switch i ~= 1
        case 1
            notConv = abs([xmCalc,ymCalc,umCalc,vmCalc]./check0 - 1) > err;
    end
    
    check0 = [xmCalc,ymCalc,umCalc,vmCalc];
    
    switch sum(sum(notConv)) == 0
        case 1
            break
    end
    
end

xOut = xmCalc;
yOut = ymCalc;
uOut = umCalc;
vOut = vmCalc;

end