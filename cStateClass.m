classdef cStateClass < handle
    
    properties
        Q
        E
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
        function obj = cStateClass(numpt, g)
            obj.Q       = zeros(3*numpt, 1);
            obj.E       = zeros(3*numpt, 1);
            
            obj.u       = zeros(numpt, 1);
            obj.rho     = zeros(numpt, 1);
            obj.p       = zeros(numpt, 1);
            obj.a       = zeros(numpt, 1);
            
            obj.rhou    = zeros(numpt, 1);
            obj.e       = zeros(numpt, 1);
            
            obj.g       = g;
        end
        
        % Get Q1 component of Q (rho)
        function Q1 = getQ1(obj)
            Q1          = obj.Q(1:3:end);
        end
        
        % Get Q2 component of Q (rho*u)
        function Q2 = getQ2(obj)
            Q2          = obj.Q(2:3:end);
        end
        
        % Get Q3 component of Q (e)
        function Q3 = getQ3(obj)
            Q3          = obj.Q(3:3:end);
        end
        
        % Get field variables (u, p, a)
        function getVar(obj)

            obj.rho     = obj.getQ1;
            obj.rhou    = obj.getQ2;
            obj.e       = obj.getQ3;

            obj.u       = obj.rhou./obj.rho;
            obj.p       = (obj.g - 1)*(obj.e - 0.5.*(obj.rhou.^2)./obj.rho);
            obj.a       = sqrt(obj.g.*obj.p./obj.rho);
            obj.H       = (obj.e + obj.p)./obj.rho;

        end
        
        % Get inviscid flux vector from Q
        function getInviscidFlux(obj)
            q1      = obj.getQ1;
            q2      = obj.getQ2;
            q3      = obj.getQ3;
            
            obj.E(1:3:end) = q2;
            obj.E(2:3:end) = (obj.g - 1)*q3 + 0.5*(3 - obj.g)*(q2.^2)./q1;
            obj.E(3:3:end) = obj.g*q3.*q2./q1 - 0.5*(obj.g - 1)*(q2.^3)./(q1.^2);
        end
    end
    
end