function [myRoots, numRoots] = rootsOfCubicEqn( a2, a1, a0 )
%ROOTSOFCUBICEQN Summary of this function goes here
%   Detailed explanation goes here
    
    g  =  a1 - (a2 * a2 / 3.0);
    h  = (2.0 * a2 * a2 * a2 - 9.0 * a2 * a1 + 27.0 * a0) / 27.0;
    gh = (g * g * g) / 27.0 + (h * h) / 4.0;
        
    if(gh > 0.0)
%      Single Real Root
        s  = -5.0e-1 * h + sqrt(gh);
        tt = -5.0e-1 * h - sqrt(gh);
            
        if (s >= 0.0)
            s  = s^(1.0/3.0);
        else
            s  = -((-s)^(1.0/3.0));
        end
        
        if (tt >= 0.0)
            tt = tt^(1.0/3.0);
        else
            tt = -((-tt)^( 1.0/3.0));    
        end
                
        myRoots  = s + tt - a2 / 3.0;
        numRoots = 1;
            
    else
%         Three Real Roots */
        xa    = -5.0e-1 * h / sqrt(-(g * g * g / 2.7e1));
        theta = (pi / 2.0 - atan(xa/(sqrt(1.0 - xa * xa))))/3.0;
            
        x1 = 2.0 * sqrt(-g / 3.0) * cos(theta);
        x2 = 2.0 * sqrt(-g / 3.0) * cos(theta + 2.0 * pi / 3.0);
        x3 = 2.0 * sqrt(-g / 3.0) * cos(theta + 4.0 * pi / 3.0);
            
        maxRoot  = max([x1,x2,x3]) - a2/3.0;

        x1mn = x1 - a2/3.0;
        x2mn = x2 - a2/3.0;
        x3mn = x3 - a2/3.0;

        if(x1mn < 0.0)
            x1mn = 1.0e10;
        end
        if(x2mn < 0.0)
            x2mn = 1.0e10;
        end
        if(x3mn < 0.0)
            x3mn = 1.0e10;
        end

        minRoot  = min([x1mn,x2mn,x3mn]);
        myRoots = [minRoot, maxRoot];
        numRoots = 3;
    end
end

