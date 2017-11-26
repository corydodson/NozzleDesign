function [output] = getprop(prop, species, varargin)
% calculates thermodynamic properties of a given species using NASA
% 9-coefficient polynomials
    
    global thermInfo
    output = 0;
    
    switch nargin
        case 3
            t = varargin{1};
        otherwise
            t = 298.15;
    end
    
    % format input
    prop = upper(prop);
    species = upper(species);
    [r,n] = size(t);
    if r < n
        r = n;
        t = t';
    end
    
    % input error check
    try 
        properties = repmat(thermInfo(species),r,1);
    catch
        disp('Species ', species,' not found in database.')
        return
    end
    
    numInt = properties{1,11};
    m = length(properties(1,:));
    greaterFilter = zeros(r, m);
    
    rUniversal = 8.3144621;
    mm = repmat(properties{1,1},r,1);
    hform = repmat(properties{1,2},r,1);
    rSpecific = rUniversal./mm;
    
    for i = 1:numInt
        propertiesFilter = double(t>=[properties{:,11 + i}]' & t<=[properties{:,12 + i}]');
        greaterFilter(:,(4 + 9*i + numInt):(4 + 9*i + numInt + 8)) = repmat(propertiesFilter,1,9);
    end
    
    properties = cell2mat(properties(:,(13 + numInt):end-2)).*greaterFilter(:,(13 + numInt):end-2);
    
    try
        properties = properties(:,1:9) + properties(:,10:18) + properties(:,19:27);
    catch
        properties = properties(:,1:9) + properties(:,10:18);
    end

%%
    
    switch prop
        case 'CP'
            t = [t.^-2, t.^-1, ones(r,1), t, t.^2, t.^3, t.^4];
            output = sum(properties(:,1:7).*t,2);
            output = output .* rSpecific;
        case 'H'
            t = [-t.^-1, log(t), t, t.^2/2, t.^3/3, t.^4/4, t.^5/5, ones(r,1)];
            output = sum(properties(:,1:8).*t,2);
            output = output .* rSpecific;
        case 'HFORM'
            output = hform/mm;
            %output = properties{2};
        case 'G'
            H = [-t.^-1, log(t), t, t.^2/2, t.^3/3, t.^4/4, t.^5/5, ones(r,1), zeros(r,1)];
            tS = [-t.^-1/2, -ones(r,1), t.*log(t), t.^2, t.^3/2, t.^4/3, t.^5/4, zeros(r,1), t.*ones(r,1)];
            t = H - tS;
            output = sum(properties(:,1:9).*t,2);
            output = output .* rSpecific;
        case 'GAMMA'
            t = [t.^-2, t.^-1, ones(r,1), t, t.^2, t.^3, t.^4];
            output = sum(properties(:,1:7).*t,2);
            output = 1./(1 - 1./output);
        case 'MM'
            output = mm;
        case 'R'
            output = rSpecific;
        case 'S'
            t = [-t.^-2/2, -t.^-1, log(t), t, t.^2/2, t.^3/3, t.^4/4, zeros(r,1), ones(r,1)];
            output = sum(properties(:,1:9).*t,2);
            output = output .* rSpecific;
        otherwise
            disp('Property not avialable.')
            return
    end

end