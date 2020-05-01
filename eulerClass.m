classdef eulerClass < handle
    
    properties
        state
        BC
        Flux
        input
    end
    
    methods
        
        % eulerStateClass Constructor
        function obj = eulerClass(input)
            obj.input           = input;
            obj.state           = cStateClass(input);
        end
        
        % Get exact boundary variables [u, rho, P, a]
        function getExactBoundaryVars(obj)
            
            Q0 = obj.BC.L.Q;
            Qe = obj.BC.R.Q;
            
            obj.BC.in.u      = Q0(2)/Q0(1);
            obj.BC.in.rho    = Q0(1);
            obj.BC.in.p      = (obj.input.thermo.g - 1)*(Q0(3) - 0.5*(Q0(2)^2)/Q0(1));
            obj.BC.in.a      = sqrt(obj.input.thermo.g*obj.BC.in.p/obj.BC.in.rho);
            
            obj.BC.out.u     = Qe(2)/Qe(1);
            obj.BC.out.rho   = Qe(1);
            obj.BC.out.p     = (obj.input.thermo.g - 1)*(Qe(3) - 0.5*(Qe(2)^2)/Qe(1));
            obj.BC.out.a     = sqrt(obj.input.thermo.g*obj.BC.out.p/obj.BC.out.rho);
            
        end
    
        % Make state non-dimensional
        function makeNonDim(obj)
            
            ndar            = [obj.input.thermo.RhoInf;                           ...
                               obj.input.thermo.RhoInf*obj.input.thermo.aInf;     ...
                               obj.input.thermo.RhoInf*obj.input.thermo.aInf^2];
                           
            obj.BC.L.Q      =  obj.BC.L.Q       ./ ndar;                   
            obj.BC.R.Q      =  obj.BC.R.Q       ./ ndar;
                                
            obj.state.Q(1:3:end)  =  obj.state.Q(1:3:end)   ./(obj.input.thermo.RhoInf);
            obj.state.Q(2:3:end)  =  obj.state.Q(2:3:end)   ./(obj.input.thermo.RhoInf*obj.input.thermo.aInf);
            obj.state.Q(3:3:end)  =  obj.state.Q(3:3:end)   ./(obj.input.thermo.RhoInf*obj.input.thermo.aInf^2);
            
        end
        
    end
    
end