function [tGuess3, pGuess3, mGuess3] = iserelimperfect(tT,pT,m,species,moleFracs)
% Calculates isentropic properties
% For a thermally perfect gas mixture with negligble chemical reaction

load thermInfo

lengthM = length(m);

switch length(tT) == lengthM
    case 0
        tT = repmat(tT,lengthM,1);
        pT = repmat(pT,lengthM,1);
end

% set up/calculate some preliminary quantities
r = mixprop('r',species,moleFracs)*1000;
ht = mixprop('h',species,moleFracs,tT);

% calculate first approximation of property values
gGuess1 = mixprop('gamma',species,moleFracs,tT);
tGuess1 = tT./(1 + (gGuess1 - 1).*m.^2./2);
hGuess1 = mixprop('h',species,moleFracs,tGuess1);
vGuess1 = sqrt(2*(ht - hGuess1)*1000);
mGuess1 = vGuess1./sqrt(gGuess1.*r.*tGuess1);

% calculate next iteration
tGuess2 = tGuess1*1.01;
gGuess2 = mixprop('gamma',species,moleFracs,tGuess2);
hGuess2 = mixprop('h',species,moleFracs,tGuess2);
vGuess2 = sqrt(2*(ht - hGuess2)*1000);
mGuess2 = vGuess2./sqrt(gGuess2.*r.*tGuess2);

% define iteration count
it = 3;
logic = abs(1 - mGuess2./m) > 10^-6;
tGuess3 = tGuess2;
gGuess3 = gGuess2;
hGuess3 = hGuess2;
vGuess3 = vGuess2;
mGuess3 = mGuess2;

% continue iteration until convergence is achieved
while sum(logic) && it<=10
    % standard property calculation block
    tGuess3(logic) = tGuess2(logic) + (tGuess2(logic) - tGuess1(logic))./(mGuess2(logic) - mGuess1(logic)).*(m(logic) - mGuess2(logic));
    gGuess3(logic) = mixprop('gamma',species,moleFracs,tGuess3(logic));
    hGuess3(logic) = mixprop('h',species,moleFracs,tGuess3(logic));
    vGuess3(logic) = sqrt(2*(ht(logic) - hGuess3(logic))*1000);
    mGuess3(logic) = vGuess3(logic)./sqrt(gGuess3(logic).*r.*tGuess3(logic));
    
    % update logic
    logic = abs(1 - mGuess3./m) > 10^-6;
   
    % shift variables and increment iteration number
    tGuess1 = tGuess2;
    mGuess1 = mGuess2;
    tGuess2 = tGuess3;
    mGuess2 = mGuess3;
    it = it + 1;
end

% p can be caluclated from s2 - s1 = -R*ln(P2/P1)
pGuess3 = pT./exp((mixprop('s',species,moleFracs,tGuess3) - mixprop('s',species,moleFracs,tT))*1000/-r);

end