classdef roeStateClass < handle
    
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Roe State Matrices
        L                       % Left State
        R                       % Right State
        Abar                    % Roe Flux Jacobian
        Lambda                  % Roe Eigenvalues
        X                       % Roe Left Eigenvectors
        Xi                      % Roe Right Eigenvectors
        flux
        Flux
        
        Abdiag
        Abdiag1
        Abdiagm1
        
        XIdiag1
        XIdiagm1
        
        Xdiag1
        Xdiagm1
        
        lam
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Thermodynamic Parameters
        g
        gm
        half_gm3
        threemg 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Numerical Parameters
        numpt
        eps
        epsi
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Block-diagonal index arrays
        dind
        oneindi
        oneindj
        twoindi
        twoindj
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sparse matrix index arrays
        Xdiag
        XIdiag
        diagp1
        diagp2
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sparse matrix constructor index arrays
        Ii
        Ij
        IAi
        IAj
    end

    methods
        
        % roeStateClass Constructor
        function obj = roeStateClass(input)
            obj.L                   = cStateClass(input.mesh.numpt+1, input.thermo.g);
            obj.R                   = cStateClass(input.mesh.numpt+1, input.thermo.g);
            
            obj.g                   = input.thermo.g;
            obj.gm                  = obj.g - 1;
            obj.half_gm3            = 0.5*(obj.g - 3);
            obj.threemg             = (3 - obj.g);
          
            obj.numpt               = input.mesh.numpt;
            obj.eps                 = input.sim.eps;
            obj.epsi                = obj.eps*ones(obj.numpt+1, 1); 
            
            obj.dind                = (1:1:(3*obj.numpt + 3))';
            obj.oneindi             = zeros(2*obj.numpt+2, 1);
            obj.oneindi(1:2:end)    = (1:3:(3*obj.numpt + 3))';
            obj.oneindi(2:2:end)    = (2:3:(3*obj.numpt + 3))';
            obj.oneindj             = obj.oneindi + ones(2*obj.numpt+2, 1);
            obj.twoindi             = (1:3:(3*obj.numpt + 3))';
            obj.twoindj             = obj.twoindi + 2*ones(obj.numpt+1, 1);  
            
            obj.Xdiag               = ones(3*obj.numpt + 3, 1);
            obj.XIdiag              = ones(3*obj.numpt + 3, 1);
            obj.diagp1              = ones(2*obj.numpt + 2, 1);
            obj.diagp2              = ones(1*obj.numpt + 1, 1);
            
            obj.Ii                  = [obj.dind;    obj.oneindi; obj.twoindi; obj.oneindj; obj.twoindj];
            obj.Ij                  = [obj.dind;    obj.oneindj; obj.twoindj; obj.oneindi; obj.twoindi];
            obj.IAi                 = [obj.oneindj; obj.oneindi; obj.oneindj; obj.twoindj];
            obj.IAj                 = [obj.oneindj; obj.oneindj; obj.oneindi; obj.twoindi];
            
            obj.lam                 = obj.Xdiag;
            obj.Xdiag1              = obj.diagp1;
            obj.Xdiagm1             = obj.diagp1;
            obj.XIdiag1             = obj.diagp1;
            obj.XIdiagm1            = obj.diagp1;
            obj.Abdiag              = obj.diagp1;
            obj.Abdiag1             = obj.diagp1;
            obj.Abdiag1(2:2:end)    = obj.gm.*obj.diagp1(2:2:end);
            obj.Abdiagm1            = obj.diagp1;
            
        end
             
        % Get Roe State
        function getState(obj)

            obj.L.getVar();
            obj.R.getVar();

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Roe Average State
            
            sqRhoL              = sqrt(obj.L.rho);
            sqRhoR              = sqrt(obj.R.rho);
            sqRhoRL             = sqRhoL + sqRhoR;
            
            rho                 = sqRhoL.*sqRhoR;
            u                   = ((obj.L.u.*sqRhoL) + (obj.R.u.*sqRhoR))./sqRhoRL;
            H                   = ((obj.L.H.*sqRhoL) + (obj.R.H.*sqRhoR))./sqRhoRL;
            a                   = sqrt(obj.gm.*(H - 0.5.*u.^2));             
            
            Lp                  = u + a;
            Lm                  = u - a;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Entropy Fix
            
            r1                  = abs(Lp) < obj.epsi;
            r2                  = abs(Lm) < obj.epsi;
            
            Lp(r1)              = 0.5.*((Lp(r1).^2)./obj.eps + obj.eps);
            Lm(r2)              = 0.5.*((Lm(r2).^2)./obj.eps + obj.eps);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate vectorized quantities

            r                   = rho       ./ (a);
            uar                 = u         .* (a);
            rHpua               = r         .* (H+uar);
            rHmua               = r         .* (H-uar);
            f                   = obj.gm    ./ (rho.*a);
            fu                  = f         .* (u);
            au2                 = a         .* (u.^2);
            u2                  = u         .^ (2);
            hu2                 = 0.5       .* (u2);
            halfr               = 0.5       .* (r);
            adgm                = a         ./ (obj.gm);
            au2dgm              = au2       ./ (obj.gm);
            Xi_11               = r         .* (adgm - hu2);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Abar

            obj.Abdiag(1:2:end)     = obj.threemg.*u;
            obj.Abdiag(2:2:end)     = obj.g.*u;
            
            obj.Abdiagm1(1:2:end)   = obj.half_gm3*u2;
            obj.Abdiagm1(2:2:end)   = H - obj.gm*u2;  
            
            obj.Abar                = sparse(obj.IAi, obj.IAj, [obj.Abdiag; obj.Abdiag1; obj.Abdiagm1; 0.5.*obj.gm.*u.^3 - u.*H]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Xi
            
            obj.XIdiag(1:3:end)     = f.*Xi_11;
            obj.XIdiag(2:3:end)     = f.*adgm - fu;
            obj.XIdiag(3:3:end)     =-f;
            
            obj.XIdiag1(1:2:end)    = fu.*r;
            obj.XIdiag1(2:2:end)    = f;
            
            obj.XIdiagm1(1:2:end)   = f.*(hu2 - au2dgm);
            obj.XIdiagm1(2:2:end)   = f.*adgm + fu;  
        
            obj.Xi                  = sparse(obj.Ii, obj.Ij, [obj.XIdiag; obj.XIdiag1; -f.*r; obj.XIdiagm1; -f.*(hu2 + au2dgm)]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % X
            
            obj.Xdiag(2:3:end)      = halfr.*Lp;
            obj.Xdiag(3:3:end)      = -0.5*rHmua;
            
            obj.Xdiag1(1:2:end)     = halfr;
            obj.Xdiag1(2:2:end)     = -halfr.*Lm;
            
            obj.Xdiagm1(1:2:end)    = u;
            obj.Xdiagm1(2:2:end)    = 0.5*rHpua;  

            obj.X                   = sparse(obj.Ii, obj.Ij, [obj.Xdiag; obj.Xdiag1; -halfr; obj.Xdiagm1; hu2]);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Lambda
        
            obj.lam(1:3:end)        = u;
            obj.lam(2:3:end)        = Lp;
            obj.lam(3:3:end)        = Lm;
            
            obj.Lambda              = sparse(obj.dind, obj.dind, obj.lam);

        end
        
        % Get Roe flux and integrated flux
        function getFlux(obj)
            
            obj.getState(); 
            obj.flux        = 0.5*(obj.Abar*(obj.L.Q + obj.R.Q) + (obj.X*abs(obj.Lambda)*obj.Xi)*(obj.L.Q - obj.R.Q)); 
            obj.Flux        = obj.flux(4:end) - obj.flux(1:end-3);
                
        end
        
    end
    
end