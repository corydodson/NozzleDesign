function [x,y,theta] = circarc(r,inp,opt)

opt = lower(opt);

switch opt
    case 'x'
        x = inp;
        x = x(x >= 0 & x <= r);
        y = 2*r - sqrt(r^2 - x.^2);
        theta = atand(x./sqrt(r^2 - x.^2));
    case 'theta'
        theta = inp - 90;
        x = r*cosd(theta);
        y = r*(2 + sind(theta));
end

theta = theta + 90;

end