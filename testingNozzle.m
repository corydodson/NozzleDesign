% boundary conditions
t0 = 3000; % total temperature [K]
p0 = 10*101325; % total pressure [Pa]
species = {'air'}; % species
moleFrac = 1; % mole fractions of species

rTd = 0.625; % radius of curvature downstream of throat
yt = 1; % half height at minimum geometric area of nozzle
d = 1; % 1 for axisymmetric, 0 for recatngular cross-section

N = 15; % number of seed points along initial value line and wall
mDesign = 5; % design Mach number
plotting = 1; % 1 for plotting on to see results in graphical form

[x,y,u,v,exitPropTable,turnPoly,strtPoly] = nozzle(t0,p0,species,moleFrac,rTd,yt,d,N,mDesign,plotting);

%%
disp(exitPropTable)

%% plot boundary geometry: initial value line and inviscid boundary
figure
plot(linspace(x(1,1),x(1,N),N),polyval(turnPoly,linspace(x(1,1),x(1,N),N)))
hold on
plot(linspace(x(1,N),x(1,end),N + 1),polyval(strtPoly,linspace(x(1,N),x(1,end),N + 1)))
axis equal
grid on
