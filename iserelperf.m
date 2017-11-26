function output = iserelperf(g,M,property)
% Finds property ratios for a given upstream mach number for isentropic flow

gP = g + 1;
gM = g - 1;

switch property
    case 't'
        % ratio of static temperatures
        output = 1 + gM.*M.^2/2;
    case 'p'
        % ratio of static pressures
        output = (1 + gM.*M.^2/2).^(g./gM);
    case 'd'
        % ratio of static desnities
        output = (1 + gM.*M.^2/2).^(1./gM);
    case 'a'
        % ratio of cross-secitonal area to critical cross-secitonal area
        output = ((2./gP).*(1 + gM.*M.^2/2)).^(gP./(2*gM))./M;
end

end