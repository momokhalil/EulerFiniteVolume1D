classdef utilitiesClass < handle
   
    properties (Constant)
        g          = 1.4;
    end
    
    methods
        
        % Get quasi-1D mach number based on S and S*
        function Mout   = getMachNumber(obj, Ar, M0)
            f = @(M) Ar - (1/M)*((2/(obj.g+1))*(1 + 0.5*(obj.g-1)*(M)^2))^((obj.g+1)/(2*(obj.g-1)));
            opts = optimoptions('fsolve', 'TolFun', 1E-14, 'TolX', 1E-14);
            options = optimset('TolX',1E-14, 'TolFun', 1E-14);
            if (Ar < 1.01)
                Mout = fsolve(f, M0, opts);
            else
                Mout = fzero(f, M0, options);
            end
        end      

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % POST PROCESSING
        
        % Make plots
        function makePlots(obj)
                obj.plotSolution();
        end
        
        % Plot numerical vs exact solutions for u and rho
        function plotSolution(obj)
            
            figure('Position', [100 10 1200 900])
            
            plot(obj.mesh.x, obj.exactSol.Ma, 'LineWidth', 1.5, 'Color', 'black');
            hold on
            plot(obj.mesh.x, obj.euler.state.u./obj.euler.state.a, '-*', 'LineWidth', 1, 'Color', 'black', 'MarkerSize', 3);

            title('Numerical and Exact solutions for Mach Number');
            xlabel('x position');
            ylabel('Mach Number');
            legend('Exact Solution', 'Numerical Solution');
            set(gca, 'FontSize', 19, 'FontName', 'CMU Sans Serif');
            pbaspect([1.3 1 1]);
            grid minor
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
            
            figure('Position', [100 10 1200 900])
            plot(obj.mesh.x, obj.exactSol.rho, 'LineWidth', 1.5, 'Color', 'black');
            hold on
            plot(obj.mesh.x, obj.euler.state.rho*obj.thermo.RhoInf, '-*', 'LineWidth', 1, 'Color', 'black', 'MarkerSize', 3);
            title('Numerical and Exact solutions for Density');
            xlabel('x position');
            ylabel('Density (kg/m^3)');
            legend('Exact Solution', 'Numerical Solution');
            set(gca, 'FontSize', 19, 'FontName', 'CMU Sans Serif');
            pbaspect([1.3 1 1]);
            grid minor
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
        end
     
        % Get exact shock tube solution
        function getExactSolution(obj)

            switch obj.sim.problemtype
                case 'shocktube'
                    aL   = sqrt(obj.thermo.g*obj.IC.pL/obj.IC.rhoL);
                    aR   = sqrt(obj.thermo.g*obj.IC.pR/obj.IC.rhoR);

                    A = sqrt(2/(obj.thermo.g*(obj.thermo.g-1)));
                    B = 2/(obj.thermo.g-1);
                    C = aL/aR;
                    alpha = (obj.thermo.g+1)/(obj.thermo.g-1);

                    f = @(P) A * ((P-1)/sqrt(1 + alpha*P)) - B*C * (1 - (obj.IC.pR*P/obj.IC.pL)^((obj.thermo.g-1)/(2*obj.thermo.g)));
                    P0 = 3;

                    Pm = fzero(f, P0);
                    
                    alpham      = (obj.thermo.g+1)/(obj.thermo.g-1);

                    p2          = Pm*obj.IC.pR;
                    p3          = p2;

                    r2          = obj.IC.rhoR*((1+alpham*Pm)/(alpham+Pm));
                    r3          = obj.IC.rhoL*(p3/obj.IC.pL)^(1/obj.thermo.g);

                    a3          = sqrt(obj.thermo.g*p3/r3);
                    a2          = sqrt(obj.thermo.g*p2/r2);

                    r           = linspace(0, obj.mesh.L, obj.mesh.numpt);
                    Mm          = linspace(0, obj.mesh.L, obj.mesh.numpt);

                    if obj.sim.nondim
                    t           = obj.sim.tfinal/obj.thermo.aInf;
                    end

                    x0          = 5;
                    V           = 2*aL*(1-(p3/obj.IC.pL)^((obj.thermo.g-1)/(2*obj.thermo.g)))/(obj.thermo.g-1);
                    C           = ((Pm-1)*aR^2)/(obj.thermo.g*V);
                    M3          = V/a3;
                    M2          = V/a2;

                    for i = 1:length(obj.mesh.x)

                        if(obj.mesh.x(i) < x0 - aL*t)

                            r(i)        = obj.IC.rhoL;
                            Mm(i)       = 0;

                        elseif(obj.mesh.x(i) > x0 - aL*t && obj.mesh.x(i) < x0 + (V*(obj.thermo.g+1)/2 - aL)*t)

                            u5          = 2*((obj.mesh.x(i) - x0)/t + aL)/(obj.thermo.g+1);
                            a5          = u5 - (obj.mesh.x(i) - x0)/t;
                            Mm(i)       = u5/a5;
                            p5          = obj.IC.pL*(a5/aL)^(2*obj.thermo.g/(obj.thermo.g-1));
                            r(i)        = obj.thermo.g*p5/a5^2;

                        elseif(obj.mesh.x(i) > x0 + (V*(obj.thermo.g+1)/2 - aL)*t && obj.mesh.x(i) < x0 + V*t)

                            r(i)        = r3;
                            Mm(i)       = M3;

                        elseif(obj.mesh.x(i) > x0 + V*t && obj.mesh.x(i) < x0 + C*t)

                            r(i)        = r2;
                            Mm(i)       = M2;

                        else

                            r(i)        = obj.IC.rhoR;
                            Mm(i)       = 0;

                        end
                    end

                    obj.exactSol.rho = r;
                    obj.exactSol.Ma = Mm;
                    
                case 'WCblast'
                    error('Error: Exact solution for problem type %s not specialized.', obj.sim.problemtype);
            end
        end
        
    end
   
end