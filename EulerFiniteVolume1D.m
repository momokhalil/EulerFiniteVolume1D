classdef EulerFiniteVolume1D < handle & utilitiesClass
    
    properties
        % General parameters
        sim, IC, input
        % Mesh
        mesh
        % Thermodynamics
        thermo
        % Euler State
        euler
        % Finite Volume Method
        fvMethod
        % Exact solution
        exactSol
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONSTRUCTOR
        
        function obj        = EulerFiniteVolume1D(input)
            
            % Input structure
            obj.input               = input;
          
            % Problem Definition
            obj.sim.problemtype     = input.sim.problemtype;
            obj.sim.ICtype          = input.sim.ICtype;
            obj.sim.BCmode          = input.sim.BCmode;
            obj.sim.BCtypeL         = input.sim.BCtypeL;
            obj.sim.BCtypeR         = input.sim.BCtypeR;
            obj.sim.nondim          = input.sim.nondim;
            obj.sim.realplot        = input.sim.realplot;
            obj.sim.timeit          = input.sim.timeit;
            obj.sim.fluxfunc        = input.sim.fluxfunc;
            obj.sim.makeplot        = input.sim.makeplot;
            obj.sim.tfinal          = input.sim.tfinal;

            % Thermodynamic Conditions
            obj.IC.pL               = input.IC.pL;
            obj.IC.pR               = input.IC.pR;
            obj.IC.rhoL             = input.IC.rhoL;
            obj.IC.rhoR             = input.IC.rhoR;
            obj.thermo.g            = input.thermo.g;
            obj.thermo.R            = input.thermo.R ;
            obj.thermo.RhoInf       = input.thermo.RhoInf;
            obj.thermo.aInf         = input.thermo.aInf;

            % Geometry
            obj.mesh.L          = input.mesh.L ;
            obj.mesh.numpt      = input.mesh.numpt;
            
            obj.createMesh();
            
            input.mesh          = obj.mesh;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Numerical Solution

            obj.sim.dtt             = input.sim.dtt;
            obj.sim.eps             = input.sim.eps;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SOLVE

            obj.sim.time.scheme     = input.sim.time.scheme;
            obj.sim.time.accurate   = input.sim.time.accurate;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Construct eulerClass
            obj.euler               = eulerClass(input);
            obj.fvMethod            = fvMethodClass(input);
            
            obj.setBCfunction();

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BOUNDARY CONDITIONS

        function setBCfunction(obj)
            switch obj.sim.BCtypeL
                case 'outflow'
                    obj.euler.BC.L.setBC = @ obj.setBCoutflowL;
                case 'reflection'
                    obj.euler.BC.L.setBC = @ obj.setBCreflectionL;
            end
            
            switch obj.sim.BCtypeR
                case 'outflow'
                    obj.euler.BC.R.setBC = @ obj.setBCoutflowR;
                case 'reflection'
                    obj.euler.BC.R.setBC = @ obj.setBCreflectionR;
            end
        end
        
        function setBoundaryConditions(obj)
            obj.euler.BC.L.setBC();
            obj.euler.BC.R.setBC();
            obj.euler.getExactBoundaryVars();              % Initialize exact boundary variables
        end

        function setBCoutflowL(obj)
            obj.euler.BC.L.Q    = obj.euler.state.Q(1:3);
        end
        
        function setBCreflectionL(obj)       
            obj.euler.BC.L.Q    = [1; -1; 1].*obj.euler.state.Q(1:3);     
        end
        
        function setBCoutflowR(obj)
            obj.euler.BC.R.Q    = obj.euler.state.Q(end-2:end);
        end
        
        function setBCreflectionR(obj)
            obj.euler.BC.R.Q    = [1; -1; 1].*obj.euler.state.Q(end-2:end);            
        end

        function updateBoundaryConditions(obj)
            obj.euler.BC.L.setBC();
            obj.euler.BC.R.setBC();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INITIALIZATION
        
        % Set initial condition
        function setInitialCondition(obj)
            
            switch obj.sim.problemtype 
                
                case 'shocktube' 
                    
                    if      strcmp(obj.sim.ICtype, 'fromBC')
                        obj.initializeStatefromBC();            % Initialize Q state from non-dimetional right boundary state
                    elseif  strcmp(obj.sim.ICtype, 'fromIC')
                        obj.initializeStatefromIC();            % Initialize Q state from non-dimetional right boundary state
                    else
                        error('IC type not specialized');
                    end
                    
                case 'WCblast' 
                    
                    if      strcmp(obj.sim.ICtype, 'fromBC')
                        obj.initializeStatefromBC();            % Initialize Q state from non-dimetional right boundary state
                    elseif  strcmp(obj.sim.ICtype, 'fromIC')
                        obj.initializeStatefromIC();            % Initialize Q state from non-dimetional right boundary state
                    else
                        error('IC type not specialized');
                    end
                    
                otherwise
                    error('Error: Problem type %s is undefined.', obj.sim.problemtype );
            end
            
            obj.euler.state.getVar();
            
        end
        
        % Make state non dimentional
        function makeNonDim(obj)
                        
            if obj.sim.nondim
                obj.euler.makeNonDim();
                obj.sim.dtt     = obj.sim.dtt   *(obj.thermo.aInf);
                obj.sim.tfinal  = obj.sim.tfinal*(obj.thermo.aInf);
            end
            
        end
        
        % Initialize interior state from specified initial condition
        function initializeStatefromIC(obj)
            
            switch obj.sim.problemtype 
                case 'shocktube' 
                    eL  = obj.IC.pL/(obj.thermo.g - 1);
                    eR  = obj.IC.pR/(obj.thermo.g - 1);

                    QL  = [obj.IC.rhoL ; 0 ; eL];
                    QR  = [obj.IC.rhoR ; 0 ; eR];
                
                    for i = 1:obj.mesh.numpt
                        if obj.mesh.x(i) < 5,   obj.euler.state.Q(3*i-2:3*i, 1) = QL;
                        else                        obj.euler.state.Q(3*i-2:3*i, 1) = QR;
                        end
                    end
                    
                case 'WCblast'
                    pL  = 1000;
                    pM  = 0.01;
                    pR  = 100;
                    
                    QL  = [1; 0; pL/(obj.thermo.g - 1)];
                    QM  = [1; 0; pM/(obj.thermo.g - 1)];
                    QR  = [1; 0; pR/(obj.thermo.g - 1)];
                    
                    for i = 1:obj.mesh.numpt
                        if obj.mesh.x(i) <= 0.1
                            obj.euler.state.Q(3*i - 2:3*i, 1) = QL;
                        elseif obj.mesh.x(i) > 0.1 && obj.mesh.x(i) <= 0.9
                            obj.euler.state.Q(3*i - 2:3*i, 1) = QM;
                        else
                            obj.euler.state.Q(3*i - 2:3*i, 1) = QR;
                        end
                    end   
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERAL FUNCTIONS        
        
        % Create geometry
        function createMesh(obj)
            
            obj.mesh.dx                 = obj.mesh.L/(obj.mesh.numpt + 1);
            obj.mesh.x                  = linspace(obj.mesh.dx, obj.mesh.L-obj.mesh.dx, obj.mesh.numpt);
            
        end
        
        % Get local time step
        function getLocalTimestep(obj)
            for i = 1:obj.mesh.numpt
                obj.sim.dt(3*i-2:3*i, 1) = [1;1;1].*obj.mesh.dx*obj.Cn/(abs(obj.euler.u(i)) + obj.euler.a(i));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TIME STEPPING FUNCTIONS
        
        % Real-time Plotting
        function realtimePlot(obj)
            figure(1)
            
            plot(obj.mesh.x, obj.euler.state.rho);
            drawnow;            
        end
        
        % Get semi-discrete form residual
        function getResidual(obj)
            obj.fvMethod.getFlux(obj.euler);
            obj.sim.Residual  = -obj.fvMethod.fluxState.Flux;
        end

        % Compute One Explicit Euler Time Step
        function explicitEulerTimeStep(obj)
            
            obj.getResidual();
            obj.euler.state.Q    = obj.euler.state.Q + obj.sim.dtt.*(obj.sim.Residual./obj.mesh.dx);   
            
        end

        % Compute One Explicit AB2 Time Step
        function AB2TimeStep(obj)
            
            obj.getResidual();
            obj.euler.state.Q    = obj.euler.state.Q + 0.5*obj.sim.dtt.*(3*obj.sim.Residual - obj.sim.Rp)./obj.mesh.dx;   
            
        end
        
        % Perform Explicit Euler First Order Time Marching
        function solveExplicitEulerTimeAccurate(obj)    
            
            obj.sim.t   = 0;    
            
            while obj.sim.t     < obj.sim.tfinal*(obj.thermo.aInf)  
                
                obj.explicitEulerTimeStep(); 
                obj.sim.t       = obj.sim.t + obj.sim.dtt;    
                
                if obj.sim.realplot
                    realtimePlot(obj);
                end
                
                obj.euler.state.getVar();    
                obj.updateBoundaryConditions();
            end    
        end  
        
        % Perform Explicit Adam's Bashforth Second Order Time Marching
        function solveExplicitAB2TimeAccurate(obj)  
            
            obj.sim.t           = 0; 
            
            obj.explicitEulerTimeStep();  
            obj.euler.state.getVar(); 
            obj.updateBoundaryConditions();
            
            obj.sim.Rp          = obj.sim.Residual;   
            obj.sim.t           = obj.sim.t + obj.sim.dtt; 
            
            while obj.sim.t     < obj.sim.tfinal
                
                obj.AB2TimeStep();   
                obj.sim.Rp      = obj.sim.Residual;  
                obj.sim.t       = obj.sim.t + obj.sim.dtt;  
                
                if obj.sim.realplot
                    realtimePlot(obj);
                end 
                
                obj.euler.state.getVar();  
                obj.updateBoundaryConditions();
            end
        end
        
        % Explicit Time Marching
        function solveExplicit(obj)
            
            obj.setInitialCondition();          % Set initial condition
            obj.setBoundaryConditions();        % Set boundary conditions
            obj.makeNonDim();                   % Non dimentionalize if specified
            obj.euler.state.getVar();         % Calculate field variables
            
            if obj.sim.timeit,  tic;            end
            
            switch obj.sim.time.accurate
                case 1
                    switch obj.sim.time.scheme
                        case 'EEuler',  obj.solveExplicitEulerTimeAccurate();
                        case 'AB2',     obj.solveExplicitAB2TimeAccurate();
                        otherwise,      error('Error: Specified time marching scheme %s is not specialized.', obj.sim.time.scheme);
                    end % switch obj.sim.time.scheme
                otherwise,              error('Error: Time accuracy mode %f is not specialized.', obj.sim.time.accurate);
            end % switch obj.sim.time.accurate
            
            if obj.sim.timeit,  toc;            end

            if obj.sim.makeplot
                obj.getExactSolution();
                obj.makePlots();
            end % makePlot
            
        end % solveExplicit
        
    end % methods
    
end % classdef


