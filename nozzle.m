function [x,y,u,v,exitPropTable,turnPoly,wallPoly] = nozzle(t0,p0,species,moleFrac,rTd,yt,d,N,mDesign,plotting)

if ~exist('thermInfo','var')
    load thermInfo
end

xDim = N + 2*(N - 1);
yDim = N + 3*(N - 1);

% calculate some thermodynamic properties
r = mixprop('r',species,moleFrac)*1000;
h0 = mixprop('h',species,moleFrac,t0);
s0 = mixprop('s',species,moleFrac,t0);

I = 10;
theta = [5*(mDesign - 1);5*mDesign;zeros(I - 2,1)];
mExit = zeros(I,1);
err = 0.01;
for i = 1:I
    
    switch i > 2
        case 1
            theta(i) = theta(i) + (theta(i) - theta(i - 1))/(mExit(i) - mExit(i - 1))*(mDesign - mExit(i));
    end
    
    [~,~,u,v,~] = kernel(t0,p0,species,moleFrac,rTd,yt,d,N,theta(i),0);
    
    vMagTemp = sqrt(u(end,end).^2 + v(end,end).^2);
    hTemp = h0 - vMagTemp.^2/2000;
    tTemp = tempfromprop(species,moleFrac,'h',hTemp);
    gTemp = mixprop('gamma',species,moleFrac,tTemp);
    mExit(i) = vMagTemp/sqrt(r*gTemp*tTemp);
    
    switch abs(mExit(i) - mDesign) < err
        case 1
            theta = theta(i);
            break
    end
    
end

[x,y,u,v,nozzleCd] = kernel(t0,p0,species,moleFrac,rTd,yt,d,N,theta,1);

turnPoly = polyfit(x(1,1:N),y(1,1:N),2);

%% calculate properties
% thermo/fluid properties
[lenRows,lenCols] = size(u);
vMag = sqrt(u.^2 + v.^2);
filt = u == 0;

h = h0 - vMag.^2/2000;
t = tempfromprop(species,moleFrac,'h',reshape(h,lenRows*lenCols,1));

g = reshape(mixprop('gamma',species,moleFrac,t),lenRows,lenCols);
g(filt) = 0;

s = reshape(mixprop('s',species,moleFrac,t),lenRows,lenCols);
s(filt) = 0;

p = p0./exp((s - s0)*1000/-r);
p(filt) = 0;

t = reshape(t,lenRows,lenCols);
t(filt) = 0;

a = sqrt(r*g.*t);
m = vMag./a;
m(filt) = 0;

% inlet mass flow rate
[tIv,pIv,~] = iserelimperfect(t0,p0,1,species,moleFrac);
gIv = mixprop('gamma',species,moleFrac,tIv);
mdotIv = pIv*sqrt(gIv/r/tIv)*(pi*yt^2)*nozzleCd;

%% fill in the rest of the nozzle
% use curve fit to find wall geometry
A = [1,-1/tand(asind(1/m(end,end)));...
     0,1];
B = [x(end,end);
    sqrt(mdotIv/pi/p(end,end)/m(end,end)*sqrt(r*t(end,end)/g(end,end)))];
X = A\B;
x(1,yDim) = X(1);
y(1,yDim) = X(2);

% match wall theta at predefined throat and theta = 0 at exit
A = [x(1,N)^3    ,x(1,N)^2   ,x(1,N)   ,1;...
    3*x(1,N)^2   ,2*x(1,N)   ,1        ,0;...
    x(1,yDim)^3  ,x(1,yDim)^2,x(1,yDim),1;...
    3*x(1,yDim)^2,2*x(1,yDim),1        ,0];
B = [y(1,N);tand(theta);y(1,yDim);tand(0)];
wallPoly = A\B;

switch plotting
    case 1
        plot(linspace(x(1,N),x(1,yDim)),polyval(wallPoly,linspace(x(1,N),x(1,yDim))),'Color',[1,0.75,0])
end

% impose straight characteristic line from last point to nozzle exit
u(:,end) = u(end,end);
v(:,end) = v(end,end);
x(:,end) = linspace(x(1,end),x(end,end),xDim);
y(:,end) = linspace(y(1,end),y(end,end),xDim);
switch plotting
    case 1
        plot(x(:,end),y(:,end),'Color',[0,0.6,0])
end

% begin finding wall nodes and subsequent nodes
for i = N:yDim - 2
    
    col = yDim;
    row = xDim - 1 - (i - N);
    
    for j = i:yDim
        
        colm1 = col - 1;
        rowm1 = row - 1;
        
        switch x(row,colm1)
            case 0
                
                % wall node
                xIn = [x(row,col);0];
                yIn = [y(row,col);0];
                uIn = [u(row,col);0];
                vIn = [v(row,col);0];
                [x(1,colm1),y(1,colm1),u(1,colm1),v(1,colm1)] = wallnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,d,wallPoly);
                
                switch plotting
                    case 1
                        plot([x(row,col),x(1,colm1)],[y(row,col),y(1,colm1)],'Color',[0,0.6,0],'MarkerEdgeColor',[0,1,0])
                end
                
                break
                
            otherwise
                
                % continue with interior nodes
                xIn = [x(row,col);x(row,colm1)];
                yIn = [y(row,col);y(row,colm1)];
                uIn = [u(row,col);u(row,colm1)];
                vIn = [v(row,col);v(row,colm1)];
                [x(rowm1,colm1),y(rowm1,colm1),u(rowm1,colm1),v(rowm1,colm1)] = internalnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,d);
                
                switch y(rowm1,colm1) > polyval(wallPoly,x(rowm1,colm1))
                    
                    case 1
                        
                        % remove old points
                        x(rowm1,colm1) = 0;
                        y(rowm1,colm1) = 0;
                        u(rowm1,colm1) = 0;
                        v(rowm1,colm1) = 0;

                        % calculate new values
                        xIn = [x(row,col);0];
                        yIn = [y(row,col);0];
                        uIn = [u(row,col);0];
                        vIn = [v(row,col);0];
                        [x(1,colm1),y(1,colm1),u(1,colm1),v(1,colm1)] = wallnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,d,wallPoly);
                        
                        switch plotting
                            case 1
                                plot([x(row,col),x(1,colm1)],[y(row,col),y(1,colm1)],'Color',[0,0.6,0],'MarkerEdgeColor',[0,1,0])
                        end
                        
                        break
                        
                end
                
                switch plotting
                    case 1
                        plot([x(row,col),x(rowm1,colm1)],[y(row,col),y(rowm1,colm1)],'Color',[0,0.6,0])
                        plot([x(row,colm1),x(rowm1,colm1)],[y(row,colm1),y(rowm1,colm1)],'Color',[0,0.6,0],'MarkerEdgeColor',[0,1,0])
                end
                
        end
        
        col = colm1;
        row = rowm1;
        
    end
end

%% calculate properties
% thermo-fluid properties
[lenRows,lenCols] = size(u);
vMag = sqrt(u.^2 + v.^2);
filt = u == 0;

h = h0 - vMag.^2/2000;
t = tempfromprop(species,moleFrac,'h',reshape(h,lenRows*lenCols,1));

g = reshape(mixprop('gamma',species,moleFrac,t),lenRows,lenCols);
g(filt) = 0;

t = reshape(t,lenRows,lenCols);
t(filt) = 0;

a = sqrt(r*g.*t);
m = vMag./a;
m(filt) = 0;

% print summary table
maxTurnAngle = theta;
velocityAtExit = vMag(end,end);
machNumberAtExit = m(end,end);
areaRatioActual = y(1,yDim)^2;
gMean = mean(g(u ~= 0));
areaRatioIsentropicPerfect = iserelperf(gMean,machNumberAtExit,'a');
exitPropTable = table(maxTurnAngle,velocityAtExit,machNumberAtExit,areaRatioActual,gMean,areaRatioIsentropicPerfect);


end