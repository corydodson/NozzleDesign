function output = tempfromprop(species,moleFracs,prop,propVal)

N = length(propVal);
tGuess = zeros(N,10);
propGuess = tGuess;

tGuess(:,1:2) = [repmat(300,N,1),repmat(2000,N,1)];
propGuess(:,1) = mixprop(prop,species,moleFracs,tGuess(:,1));
propGuess(:,2) = mixprop(prop,species,moleFracs,tGuess(:,2));

it = 3;
logic = abs(propGuess(:,it)./propVal - 1) > 1e-6;

while sum(logic) && it<=10
    
    tGuess(logic,it) = tGuess(logic,it - 1) + (tGuess(logic,it - 1) - tGuess(logic,it - 2)) ./ (propGuess(logic,it - 1) - propGuess(logic,it - 2)) .* (propVal(logic) - propGuess(logic,it - 1));
    tGuess(~logic,it) = tGuess(~logic,it - 1);
    propGuess(logic,it) = mixprop(prop,species,moleFracs,tGuess(logic,it));
    propGuess(~logic,it) = propGuess(~logic,it - 1);
    
    logic = abs(propGuess(:,it)./propVal - 1) > 1e-6;
    it = it + 1;
    
end

output = tGuess(:,it - 1);

end