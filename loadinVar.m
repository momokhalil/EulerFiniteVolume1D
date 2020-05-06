function input = loadinVar()

    % Problem Definition
    input.sim.problemtype   = 'shocktube';
    input.sim.ICtype        = 'fromIC';
    input.sim.BCtypeL       = 'outflow';
    input.sim.BCtypeR       = 'outflow';
    input.sim.BCmode        = 'dirichlet';
    input.sim.fluxfunc      = 'Roe';
    input.sim.realplot      = 1;
    input.sim.nondim        = 1;
    input.sim.makeplot      = 1;
    input.sim.timeit        = 0;
    input.sim.tfinal        = 6.1/1000;
    input.sim.dtt           = 0.0025/1000;
    input.sim.eps           = 0.000001;
    input.sim.time.scheme   = 'RK2';
    input.sim.time.accurate = 1;
    input.order             = 2;
    
    % Initial Condition
    input.IC.pL             = 10^5;
    input.IC.pR             = 10^4;
    input.IC.rhoL           = 1;
    input.IC.rhoR           = 0.125;
    
    % Thermodynamic Conditions
    input.thermo.g          = 1.4;
    input.thermo.R          = 287;
    input.thermo.RhoInf     = 1;
    input.thermo.aInf       = 343;

    % Geometry
    input.mesh.L            = 10;
    input.mesh.numpt        = 500;
    
    % Limiter type
    input.limiter.type      = 'vanAlbada';

end