% boundary conditions
t0 = 3000;
p0 = 10*101325;
species = {'air'};
moleFrac = 1;

rTd = 0.625;
yt = 1;
d = 1;

N = 15;
mDesign = 5;
plotting = 1;

tic
[x,y,u,v,exitPropTable,turnPoly,strtPoly] = nozzle(t0,p0,species,moleFrac,rTd,yt,d,N,mDesign,plotting);
toc

%%
disp(exitPropTable)

%%
figure
plot(linspace(x(1,1),x(1,N),N),polyval(turnPoly,linspace(x(1,1),x(1,N),N)))
hold on
plot(linspace(x(1,N),x(1,end),N + 1),polyval(strtPoly,linspace(x(1,N),x(1,end),N + 1)))
axis equal
grid on