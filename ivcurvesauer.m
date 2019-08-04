function [x,y,u,v] = ivcurvesauer(t0,p0,species,moleFracs,d,rt,yt,N)

r = mixprop('r',species,moleFracs)*1000;
h0 = mixprop('h',species,moleFracs,t0);

[t, ~, ~] = iserelimperfect(t0,p0,repmat(1.2,N,1),species,moleFracs);
g = mixprop('gamma',species,moleFracs,t);
vMag = sqrt(r.*g.*t);

I = 20;
err = 1e-4;
for i = 1:I
    alpha = sqrt((1 + d)./(g + 1)/rt/yt);
    ep = -yt/(2*(3 + d))*sqrt((g + 1)*(1 + d)/(rt/yt));

    y = linspace(yt,0,N)';
    x = -(g + 1).*alpha.*y.^2/2/(3 + d);

    [u,v] = velsaur(x,y,d,alpha,g);
    a = sqrt(r*g.*t);
    u = (1 + u).*a;
    v = v.*a;
    vTemp = sqrt(u.^2 + v.^2);
    
    notConv = abs(vTemp./vMag - 1) > err;
    switch sum(notConv) == 0
        case 1
            break
    end
    
    vMag = vTemp;
    h = h0 - vMag.^2/2000;
    t = tempfromprop(species,moleFracs,'h',h);
    g = mixprop('gamma',species,moleFracs,t);

end

x = x - ep(end);

end