function [xOut,yOut,uOut,vOut] = internalnode(t0,species,moleFracs,x,y,u,v,d)

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

sp = vm;

N = 20;
err = 1e-4;
notConv = true(1,4);
for i = 1:N
    
    ymCalc = (ym + ymCalc)/2;
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
    ap = a(plus);
    
    mu = asind(a./vMag);
    mum = mu(minus);
    mup = mu(plus);
    
    theta = atand(vCalc./uCalc);
    thetam = theta(minus);
    thetap = theta(plus);
    
    lm = tand(thetam - mum);
    lp = tand(thetap + mup);
    
    q = uCalc.^2 - a.^2;
    qm = q(minus);
    qp = q(plus);
    
    rm = 2*umCalc.*vmCalc - qm.*lm;
    rp = 2*upCalc.*vpCalc - qp.*lp;
    
    sm = d*am.^2.*vmCalc./ymCalc;
    switch isempty(sp(ypCalc == 0))
        case 1
            sp = d*ap.^2.*vpCalc./ypCalc;
        otherwise
            sp(ypCalc ~= 0) = d*ap(ypCalc ~= 0).^2.*vpCalc(ypCalc ~= 0)./ypCalc(ypCalc ~= 0);
            sp(ypCalc == 0) = sm(end);
    end
    
    A = zeros(L);
    B = zeros(L,1);
    for j = 1:L - 1
        A(2*(j - 1) + 1:2*j,2*(j - 1) + 1:2*j) = [1,-lp(j);...
                                                  1,-lm(j)];
        B(2*(j - 1) + 1:2*j) = [yp(j) - lp(j)*xp(j);...
                                ym(j) - lm(j)*xm(j)];
    end
    
    X = A\B;
    ymCalc = X(1:2:end);
    ypCalc = ymCalc;
    xmCalc = X(2:2:end);
    xpCalc = xmCalc;
    
    for j = 1:L - 1
        A(2*(j - 1) + 1:2*j,2*(j - 1) + 1:2*j) = [qp(j),rp(j);...
                                                  qm(j),rm(j)];
        B(2*(j - 1) + 1:2*j) = [sp(j)*(xpCalc(j) - xp(j)) + qp(j)*up(j) + rp(j)*vp(j);...
                                sm(j)*(xmCalc(j) - xm(j)) + qm(j)*um(j) + rm(j)*vm(j)];
    end
    
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
    
    switch isnan(xmCalc) | isnan(ymCalc) | isnan(umCalc) | isnan(vmCalc)
        case 1
            break
    end
    
end

xOut = xmCalc;
yOut = ymCalc;
uOut = umCalc;
vOut = vmCalc;

end