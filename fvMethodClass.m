classdef fvMethodClass < handle
    
    properties
        order
        fluxState
        limiter
        input
    end
    
    methods
        
        % fvMethodClass Constructor
        function obj    = fvMethodClass(input)
            obj.input   = input;
            obj.order   = input.order;
            
            obj.limiter.type = obj.input.limiter.type;
            
            obj.setFluxFunction();
            obj.setFluxLimiterFunction();
            obj.setLimitedStatefunction();
            obj.setSlopeCorrectionFunction();
        end
        
        % Select flux function type
        function setFluxFunction(obj)
            switch obj.input.sim.fluxfunc
                case 'Roe'
                    obj.fluxState = roeStateClass(obj.input);
            end
        end
        
        % Set flux limiter function
        function setFluxLimiterFunction(obj) 
            switch obj.input.limiter.type
                case 'vanAlbada',   obj.limiter.func    = @(x) (x.^2 + x)./(x.^2 + 1);
                case 'vanLeer',     obj.limiter.func    = @(x) (x + abs(x))./(x + 1);
            end
        end
        
        % Select function to either limit solution or leave unlimited
        function setLimitedStatefunction(obj)
            switch obj.limiter.type
                case 'none'
                    obj.limiter.getState  = @ obj.reconstructUnlimitedState;
                otherwise
                    obj.limiter.getState  = @ obj.reconstructLimitedState;
            end
        end
        
        % Get slope used to calculate limiter function
        function getSlope(obj, eulerState)
            
            Q       = eulerState.state.Q;
            
            rrho    = zeros(obj.input.mesh.numpt, 1);
            rrhou   = zeros(obj.input.mesh.numpt, 1);
            re      = zeros(obj.input.mesh.numpt, 1);
            
            if      Q(1) - eulerState.BC.L.Q(1) == 0 && Q(4) - Q(1) > 0,            rrho(1) = 1000000;
            elseif  Q(1) - eulerState.BC.L.Q(1) == 0 && Q(4) - Q(1) <= 0,           rrho(1) = 0;
            else    rrho(1)     = (Q(4) - Q(1))/(Q(1) - eulerState.BC.L.Q(1));
            end

            if      Q(2) - eulerState.BC.L.Q(2) == 0 && Q(5) - Q(2) > 0,            rrhou(1) = 1000000;
            elseif  Q(2) - eulerState.BC.L.Q(2) == 0 && Q(5) - Q(2) <= 0,           rrhou(1) = 0;
            else    rrhou(1)    = (Q(5) - Q(2))/(Q(2) - eulerState.BC.L.Q(2));
            end   

            if      Q(3) - eulerState.BC.L.Q(3) == 0 && Q(6) - Q(3) > 0,            re(1) = 1000000;
            elseif  Q(3) - eulerState.BC.L.Q(3) == 0 && Q(6) - Q(3) <= 0,           re(1) = 0;
            else    re(1)       = (Q(6) - Q(3))/(Q(3) - eulerState.BC.L.Q(3));
            end  
                
            for i = 2:obj.input.mesh.numpt-1
                
                if      (Q(3*i-2) - Q(3*i-5) == 0) && (Q(3*i+1) - Q(3*i-2) > 0) ,    rrho(i) = 1000000;
                elseif  (Q(3*i-2) - Q(3*i-5) == 0) && (Q(3*i+1) - Q(3*i-2) <= 0) ,   rrho(i) = 0; 
                else     rrho(i) = (Q(3*i+1) - Q(3*i-2))/(Q(3*i-2) - Q(3*i-5));
                end
                
                if      (Q(3*i-1) - Q(3*i-4) == 0) && (Q(3*i+2) - Q(3*i-1) > 0) ,    rrhou(i) = 1000000;
                elseif  (Q(3*i-1) - Q(3*i-4) == 0) && (Q(3*i+2) - Q(3*i-1) <= 0) ,   rrhou(i) = 0;
                else     rrhou(i) = (Q(3*i+2) - Q(3*i-1))/(Q(3*i-1) - Q(3*i-4));
                end   
                
                if      (Q(3*i+0) - Q(3*i-3) == 0) && (Q(3*i+3) - Q(3*i+0) > 0) ,    re(i) = 1000000;
                elseif  (Q(3*i+0) - Q(3*i-3) == 0) && (Q(3*i+3) - Q(3*i+0) <= 0) ,   re(i) = 0;
                else     re(i) = (Q(3*i+3) - Q(3*i+0))/(Q(3*i+0) - Q(3*i-3));
                end                   
            end
           
            if      Q(end-2) - Q(end-5) == 0 && (eulerState.BC.R.Q(1) - Q(end-2)) > 0,  rrho(end) = 100000;
            elseif  Q(end-2) - Q(end-5) == 0 && (eulerState.BC.R.Q(1) - Q(end-2)) <= 0, rrho(end) = 0;
            else    rrho(end)     = (eulerState.BC.R.Q(1) - Q(end-2))/(Q(end-2) - Q(end-5));
            end

            if      Q(end-1) - Q(end-4) == 0 && (eulerState.BC.R.Q(2) - Q(end-1)) > 0,  rrhou(end) = 100000;
            elseif  Q(end-1) - Q(end-4) == 0 && (eulerState.BC.R.Q(2) - Q(end-1)) <= 0, rrhou(end) = 0;
            else    rrhou(end)     = (eulerState.BC.R.Q(2) - Q(end-1))/(Q(end-1) - Q(end-4));
            end 

            if      Q(end-0) - Q(end-3) == 0 && (eulerState.BC.R.Q(3) - Q(end-0)) > 0,  re(end) = 100000;
            elseif  Q(end-0) - Q(end-3) == 0 && (eulerState.BC.R.Q(3) - Q(end-0)) <= 0, re(end) = 0;
            else    re(end)       = (eulerState.BC.R.Q(3) - Q(end))/(Q(end) - Q(end-3));
            end  
                    
            rrho(rrho == -1)     = -1 + obj.input.sim.eps;
            rrhou(rrhou == -1)   = -1 + obj.input.sim.eps;
            re(re == -1)         = -1 + obj.input.sim.eps;
            
            obj.limiter.slope.rho   = rrho;
            obj.limiter.slope.rhou  = rrhou;
            obj.limiter.slope.e     = re;            
        end
                
        % Get flux limiter
        function getFluxLimiter(obj)
            
            obj.limiter.Prho            = zeros(obj.input.mesh.numpt, 1);
            obj.limiter.Prhou           = zeros(obj.input.mesh.numpt, 1);
            obj.limiter.Pe              = zeros(obj.input.mesh.numpt, 1); 
            
            obj.limiter.Prho    (obj.limiter.slope.rho > 0)     = obj.limiter.func(obj.limiter.slope.rho(obj.limiter.slope.rho > 0));            
            obj.limiter.Prhou   (obj.limiter.slope.rhou > 0)    = obj.limiter.func(obj.limiter.slope.rhou(obj.limiter.slope.rhou > 0));               
            obj.limiter.Pe      (obj.limiter.slope.e > 0)       = obj.limiter.func(obj.limiter.slope.e(obj.limiter.slope.e > 0));    

            obj.limiter.slope.correction();

            obj.limiter.Phi             = zeros(3*obj.input.mesh.numpt, 1);
            obj.limiter.Phi(1:3:end)    = obj.limiter.Prho;
            obj.limiter.Phi(2:3:end)    = obj.limiter.Prhou;
            obj.limiter.Phi(3:3:end)    = obj.limiter.Pe;   
            
        end
        
        function setSlopeCorrectionFunction(obj)
            switch obj.order
                case 1, obj.limiter.slope.correction    = @ obj.correctFirstOrder;
                case 2, obj.limiter.slope.correction    = @ obj.correctSecondOrder;
            end
        end
        
        function correctFirstOrder(obj)
            obj.limiter.Prho            = obj.limiter.Prho;
            obj.limiter.Prhou           = obj.limiter.Prhou;
            obj.limiter.Pe              = obj.limiter.Pe; 
        end
        
        function correctSecondOrder(obj)
            obj.limiter.Prho            = (0.5./(obj.limiter.slope.rho + 1)).*obj.limiter.Prho;
            obj.limiter.Prhou           = (0.5./(obj.limiter.slope.rhou + 1)).*obj.limiter.Prhou;
            obj.limiter.Pe              = (0.5./(obj.limiter.slope.e + 1)).*obj.limiter.Pe; 
        end
        
        % Get unlimited Left and Right Interface States
        function reconstructUnlimitedState(obj, eulerState)
            Qq                  = [eulerState.BC.L.Q; eulerState.state.Q; eulerState.BC.R.Q]; 
            PQ                  = (Qq(7:end) - Qq(1:end-6));
            
            if obj.order == 1
                obj.fluxState.L.Q   = Qq(1:end-3);    
                obj.fluxState.R.Q   = Qq(4:end);    
            else
                obj.fluxState.L.Q   = Qq(1:end-3) + [[0;0;0]; PQ];    
                obj.fluxState.R.Q   = Qq(4:end)   - [PQ; [0;0;0]];                
            end
            
        end
        
        % Get limited Left and Right Interface States
        function reconstructLimitedState(obj, eulerState)
            obj.getSlope(eulerState);            
            obj.getFluxLimiter();  

            Qq                  = [eulerState.BC.L.Q; eulerState.state.Q; eulerState.BC.R.Q];          
            PQ                  = obj.limiter.Phi.*(Qq(7:end) - Qq(1:end-6));

            obj.fluxState.L.Q   = Qq(1:end-3) + (obj.order-1)*[[0;0;0]; PQ];    
            obj.fluxState.R.Q   = Qq(4:end)   - (obj.order-1)*[PQ; [0;0;0]];
            
        end
        
        % Get flux
        function getFlux(obj, eulerState)
            
            obj.limiter.getState(eulerState);
            obj.fluxState.getFlux();
            
        end

    end

end