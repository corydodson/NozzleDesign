function [output] = mixProp(prop, species, moleFractions, varargin)

% format input
prop = upper(prop);
species = upper(species);
[m,n] = size(moleFractions);
if n > m
    moleFractions = moleFractions';
end
moleFractions = moleFractions./repmat(sum(moleFractions,1),m,1);

% set initial values
output = 0;
mixMolarMass = 0;
rUniversal = 8.3144621;

% error check: same length input arrays
switch length(species)
    case length(moleFractions)
    otherwise
        disp('Species and mole fraction arrays not of equal length.')
        return
end

% error check: set default temperature if none present
switch nargin
    case 4
        t = varargin{1};
    otherwise
        t = ones(m,n)*298.15;
end

% error check: make sure property is available
switch prop
    case 'CP'
    case 'G'
    case 'H'
    case 'HFORM'
    case 'GAMMA'
    case 'MM'
        for i = 1:length(species)
            output = output + getprop('MM',species{i}) * moleFractions(i);
        end
        return
    case 'R'
        for i = 1:length(species)
            output = output + getprop('MM',species{i}) * moleFractions(i);
        end
        output = rUniversal/output;
        return
    case 'S'
    otherwise
        disp('Property not avialable.')
        return
end

% use ideal gas mixing rule to calculate properties

for i = 1:length(species)
    mmTemp = getprop('MM',species{i});
    mixMolarMass = mixMolarMass + mmTemp * moleFractions(i);
    output = output + getprop(prop,species{i},t) * mmTemp * moleFractions(i);
end

% convert from mole to mass based output
output = output / mixMolarMass;

end