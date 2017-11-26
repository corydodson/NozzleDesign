function [x,y,u,v,nozzleCd] = kernel(t0,p0,species,moleFrac,rTd,yt,d,N,thetaMax,plotting)

% pre-allocate arrays
xDim = N + 2*(N - 1);
yDim = N + 3*(N - 1);
x = zeros(xDim,yDim);
y = x;
u = x;
v = x;
t = x;
g = x;
a = x;

%% Calculate initial value line
% define initial value line
[t(1:N,1),~,mxIv,myIv,~,~,x(1:N,1),y(1:N,1),nozzleCd] = ivcurvekliegl(t0,p0,species,moleFrac,rTd,yt,N);
g(1:N,1) = mixprop('gamma',species,moleFrac,t(1:N,1));
r = mixprop('r',species,moleFrac)*1000;
a(1:N,1) = sqrt(r*g(1:N,1).*t(1:N,1));
u(1:N,1) = mxIv.*a(1:N,1);
v(1:N,1) = myIv.*a(1:N,1);

switch plotting
    case 1
        figure
        plot(x(1:N,1),y(1:N,1),'Color',[0,0.75,1])

        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + 1.5*ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - 1.5*(ti(1) + ti(3));
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        hold on
        grid on
end

%% define wall just downstream of throat
[xTd,yTd,thetaTd] = circarc(rTd,linspace(0,thetaMax,N),'theta');

x(1,1:N) = xTd/rTd;
y(1,1:N) = yTd/rTd;

switch plotting
    case 1
        plot(x(1,1:N),y(1,1:N),'Color',[1,0.75,0])
        axis equal
end

%% calculate nodes from initial value (IV) line
count = 1;
for i = 1:N - 1
    
    ind = N - i;
    
    for j = 1:count
        
        jp1 = j + 1;
        jp2 = j + 2;
        
        ind = ind + 1;
        indP1 = ind + 1;
        indM1 = ind - 1;
        
        % calculate interior points
        xIn = x(indM1:ind,j);
        yIn = y(indM1:ind,j);
        uIn = u(indM1:ind,j);
        vIn = v(indM1:ind,j);
        [x(ind,jp1),y(ind,jp1),u(ind,jp1),v(ind,jp1)] = internalnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,d);
        
        switch plotting
            case 1
                plot([x(indM1,j),x(ind,jp1)],[y(indM1,j),y(ind,jp1)],'Color',[0,0.6,0])
                plot([x(ind,j),x(ind,jp1)],[y(ind,j),y(ind,jp1)],'Color',[0,0.6,0],'MarkerEdgeColor',[0,1,0])
        end
        
    end
    
    % calculate axis node
    xIn = [x(ind,jp1);x(ind,jp1)];
    yIn = [y(ind,jp1);-y(ind,jp1)];
    uIn = [u(ind,jp1);u(ind,jp1)];
    vIn = [v(ind,jp1);-v(ind,jp1)];
    [x(indP1,jp2),y(indP1,jp2),u(indP1,jp2),v(indP1,jp2)] = internalnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,d);
    
    switch plotting
        case 1
            plot([x(ind,jp1),x(indP1,jp2)],[y(ind,jp1),y(indP1,jp2)],'Color',[0,0.6,0],'MarkerEdgeColor',[0,1,0])
    end
    
    count = count + 2;
    
end

%% calculate nodes from wall remaining nodes in kernel
minRow = 2;
skips = 0;
for i = 1:N - 1
    
    ip1 = i + 1;
    ip2 = i + 2;
    
    % wall node
    xIn = [x(1,i);x(minRow,ip1 + skips);0;x(1,ip1)];
    yIn = [y(1,i);y(minRow,ip1 + skips);0;y(1,ip1)];
    uIn = [u(1,i);u(minRow,ip1 + skips);0;0];
    vIn = [v(1,i);v(minRow,ip1 + skips);0;0];
    thetaIn = [0;0;0;thetaTd(ip1)];
    try
        [x(1,ip1),y(1,ip1),u(1,ip1),v(1,ip1)] = invwallnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,thetaIn,d);
    catch
        xIn = [x(1,i);x(minRow + 1,ip1 + skips + 1);0;x(1,ip1)];
        yIn = [y(1,i);y(minRow + 1,ip1 + skips + 1);0;y(1,ip1)];
        uIn = [u(1,i);u(minRow + 1,ip1 + skips + 1);0;0];
        vIn = [v(1,i);v(minRow + 1,ip1 + skips + 1);0;0];
        [x(1,ip1),y(1,ip1),u(1,ip1),v(1,ip1)] = invwallnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,thetaIn,d);
    end
    
    % near wall node
    success = 0;
    while success == 0
        
        xIn = [x(1,ip1);x(minRow,ip1 + skips)];
        yIn = [y(1,ip1);y(minRow,ip1 + skips)];
        uIn = [u(1,ip1);u(minRow,ip1 + skips)];
        vIn = [v(1,ip1);v(minRow,ip1 + skips)];
        [x(minRow,ip2 + skips),y(minRow,ip2 + skips),u(minRow,ip2 + skips),v(minRow,ip2 + skips)] = internalnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,d);
        
        % discontinue lines that would reflect in initial expansion region
        filt = x(minRow,ip2 + skips) < x(1,N) & (imag(x(minRow,ip2 + skips)) ~= 0 | y(minRow,ip2 + skips) > interp1(x(1,1:N),y(1,1:N),real(x(minRow,ip2 + skips))));
        switch sum(filt) >= 1
            case 1
                x(minRow,ip2 + skips) = 0;
                y(minRow,ip2 + skips) = 0;
                u(minRow,ip2 + skips) = 0;
                v(minRow,ip2 + skips) = 0;
                minRow = minRow + 1;
                skips = skips + 1;
            otherwise
                success = 1;
                
                switch plotting
                    case 1
                        plot([x(1,ip1),x(minRow,ip2 + skips)],[y(1,ip1),y(minRow,ip2 + skips)],'Color',[0,0.6,0])
                        plot([x(minRow,ip1 + skips),x(minRow,ip2 + skips)],[y(minRow,ip1 + skips),y(minRow,ip2 + skips)],'Color',[0,0.6,0],'MarkerEdgeColor',[0,1,0])
                end
                
        end
    end

    % continue with interior nodes
    ind = i + 2 + skips;
    for j = minRow:count - 1
        
        jp1 = j + 1;
        jp2 = j + 2;
        indP1 = ind + 1;
        indP2 = ind + 2;

        xIn = [x(j,ind);x(jp1,ind)];
        yIn = [y(j,ind);y(jp1,ind)];
        uIn = [u(j,ind);u(jp1,ind)];
        vIn = [v(j,ind);v(jp1,ind)];
        [x(jp1,indP1),y(jp1,indP1),u(jp1,indP1),v(jp1,indP1)] = internalnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,d);
        
        switch plotting
            case 1
                plot([x(j,ind),x(jp1,indP1)],[y(j,ind),y(jp1,indP1)],'Color',[0,0.6,0])
                plot([x(jp1,ind),x(jp1,indP1)],[y(jp1,ind),y(jp1,indP1)],'Color',[0,0.6,0],'MarkerEdgeColor',[0,1,0])
        end

        ind = ind + 1;
    end

    % calculate axis node
    xIn = [x(jp1,indP1);x(jp1,indP1)];
    yIn = [y(jp1,indP1);-y(jp1,indP1)];
    uIn = [u(jp1,indP1);u(jp1,indP1)];
    vIn = [v(jp1,indP1);-v(jp1,indP1)];
    [x(jp2,indP2),y(jp2,indP2),u(jp2,indP2),v(jp2,indP2)] = internalnode(t0,species,moleFrac,xIn,yIn,uIn,vIn,d);
    
    switch plotting
        case 1
            plot([x(jp1,indP1),x(jp2,indP2)],[y(jp1,indP1),y(jp2,indP2)],'Color',[0,0.6,0],'MarkerEdgeColor',[0,1,0])
    end

    count = count + 1;
    
end

end