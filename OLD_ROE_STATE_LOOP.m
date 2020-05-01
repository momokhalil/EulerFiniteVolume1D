            for i = 1:obj.numpt +1    
                obj.Abar(3*i-2:3*i, 3*i-2:3*i)     = [0,                                1,                      0;
                                                      obj.half_gm3*u2(i),               obj.threemg*u(i),       obj.gm;
                                                      0.5*obj.gm*u(i)^3 - u(i)*Hr(i),   Hr(i) - obj.gm*u2(i),   obj.g*u(i)];

                obj.Lambda(3*i-2, 3*i-2)           =  u(i);
                obj.Lambda(3*i-1, 3*i-1)           =  Lp(i);
                obj.Lambda(3*i-0, 3*i-0)           =  Lm(i);

                obj.X(3*i-2:3*i, 3*i-2:3*i)        = [1,                                halfr(i),              -halfr(i);
                                                      u(i),                             halfr(i)*Lp(i),        -halfr(i)*Lm(i);
                                                      0.5*u2(i),                        0.5*rHpua(i),          -0.5*rHmua(i)];

                obj.Xi(3*i-2:3*i, 3*i-2:3*i)  = f(i)*[r(i)*(ar(i)/obj.gm-0.5*u2(i)),    r(i)*u(i),             -r(i);
                                                      0.5*u2(i) - au2(i)/obj.gm,        ar(i)/obj.gm - u(i),    1;
                                                     -0.5*u2(i) - au2(i)/obj.gm,        ar(i)/obj.gm + u(i),   -1];

            end  