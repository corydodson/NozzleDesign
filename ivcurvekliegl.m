function [t,p,u,v,x,y,z,r,cd] = ivcurvekliegl(t0,p0,species,moleFracs,rt,yt,N)
%%
% inputs
%     t0: total (stagnation) temperature [K]
%     p0: total (stagnation) pressure [Pa]
%     species: vertical vector (1-D) cell array of strings with gaseous species names
%     moleFracs: vertical vector (1-D) array of doubles with corresponding gaseous species mole fractions
%     rt: radius of curvature of the wall at the axial position of the geometric throat
%     yt: distance between the nozzle centerline and the wall at the axial position of the geometric throat
%     N: number of points along the initial value line to calculate properties
% outputs
%     t: static temperature at points
%     p: static pressure at points
%     x: x-location of points
%     y: y-location of points
%     u: Mach number of the axial velocity component
%     v: Mach number of the radial velocity component

% preallocate variables
I = 10;
zGuess = zeros(N,I);
uGuess = zGuess;
vGuess = uGuess;
mGuess = vGuess;
r = linspace(1,0,N)';
r = [r(1:end - 1);-r(end - 1)];

% calculate thermodynamic properties
rSpec = mixprop('r',species,moleFracs)*1000;
[tCrit, ~, ~] = iserelimperfect(t0,p0,1,species,moleFracs);
gCrit = mixprop('gamma',species,moleFracs,tCrit);

Rt = rt/yt;
cd = 1 - (gCrit + 1)/(1 + Rt)^2.*(1/96 - (8*gCrit - 27)/2304/(1 + Rt) + (754*gCrit^2 - 757*gCrit + 3633)/276480/(1 + Rt)^2);

% begin guessing
i = 1;
zGuess(:,i) = 0;
[uGuess(:,i),vGuess(:,i),~,~] = velkliegl(zGuess(:,i),r,rt,yt,gCrit);
mGuess(:,i) = sqrt(uGuess(:,i).^2 + vGuess(:,i).^2);
[tCurve,~,~] = iserelimperfect(t0,p0,mGuess(:,i),species,moleFracs);
gCurve = mixprop('gamma',species,moleFracs,tCurve);
% plot(zGuess(:,i),r);hold on

i = 2;
zGuess(:,i) = 0.2;
[uGuess(:,i),vGuess(:,i),~,~] = velkliegl(zGuess(:,i),r,rt,yt,gCurve);
mGuess(:,i) = sqrt(uGuess(:,i).^2 + vGuess(:,i).^2);
[tCurve,~,~] = iserelimperfect(t0,p0,mGuess(:,i),species,moleFracs);
gCurve = mixprop('gamma',species,moleFracs,tCurve);
% plot(zGuess(:,i),r);

notconv = abs(vGuess(:,i)) > 1e-6;
for i = 3:I
    
    zGuess(notconv,i) = zGuess(notconv,i - 1) + (zGuess(notconv,i - 1) - zGuess(notconv,i - 2))./(vGuess(notconv,i - 1) - vGuess(notconv,i - 2)).*(0 - vGuess(notconv,i - 1));
    zGuess(~notconv,i) = zGuess(~notconv,i - 1);
    [uGuess(:,i),vGuess(:,i),x,y] = velkliegl(zGuess(:,i),r,rt,yt,gCurve);
    mGuess(:,i) = sqrt(uGuess(:,i).^2 + vGuess(:,i).^2);
    [tCurve,~,~] = iserelimperfect(t0,p0,mGuess(:,i),species,moleFracs);
    gCurve = mixprop('gamma',species,moleFracs,tCurve);
    
    notconv = abs(vGuess(:,i)) > 1e-6 & vGuess(:,i) > 0;
    notconv(end) = notconv(end - 1);
%     plot(zGuess(:,i),r);
    switch sum(notconv)
        case 0
            break
    end
    
end

u = uGuess(:,i);
v = vGuess(:,i);
z = zGuess(:,i);
t = tCurve;

z(end) = interp1(r,z,0,'poly2');
u(end) = griddata(r,z,u,0,z(end));
v(end) = griddata(r,z,v,0,z(end));
t(end) = griddata(r,z,t,0,z(end));
r(end) = 0;

s0 = mixprop('s',species,moleFracs,t0);
s = mixprop('s',species,moleFracs,tCurve);
p = p0./exp((s - s0)/-rSpec);

end