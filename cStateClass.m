classdef cStateClass < handle
    
    properties
        Q
        rho
        rhou
        e
        u
        p
        a
        H
        
        g
    end
    
    methods
        % cStateClass constructor
        function obj = cStateClass(input)
            obj.Q       = zeros(3*input.mesh.numpt, 1);
            obj.u       = zeros(input.mesh.numpt, 1);
            obj.rho     = zeros(input.mesh.numpt, 1);
            obj.p       = zeros(input.mesh.numpt, 1);
            obj.a       = zeros(input.mesh.numpt, 1);
            
            obj.g       = input.thermo.g;
        end
        
        % Get Q1
        function Q1 = getQ1(obj)
            Q1          = obj.Q(1:3:end);
        end
        
        % Get Q2
        function Q2 = getQ2(obj)
            Q2          = obj.Q(2:3:end);
        end
        
        % Get Q3
        function Q3 = getQ3(obj)
            Q3          = obj.Q(3:3:end);
        end
        
        % Get field variables
        function getVar(obj)

            obj.rho     = obj.getQ1;
            obj.rhou    = obj.getQ2;
            obj.e       = obj.getQ3;

            obj.u       = obj.rhou./obj.rho;
            obj.p       = (obj.g - 1)*(obj.e - 0.5.*(obj.rhou.^2)./obj.rho);
            obj.a       = sqrt(obj.g.*obj.p./obj.rho);
            obj.H       = (obj.e + obj.p)./obj.rho;

        end
        
    end
    
end