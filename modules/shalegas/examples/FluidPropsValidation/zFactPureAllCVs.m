% '------------------------------------------------------------------------
% 'Calculation of Z Factor for Pure substance with PR-EOS
% '------------------------------------------------------------------------
function z = zFactPureAllCVs(pres,TinK,pCrit,tCrit,accFact,RR,numCVs)

    coder.inline('never');
    z = zeros(numCVs,1);
    
    for i=1:numCVs
        Tr     = TinK/tCrit;
        %applies to accFact from 0 to 2. Robinson et al.,1985.
        if (accFact<0.5)
            kappa  = (-0.26992*accFact + 1.54226)*accFact + 0.37464;
        else
            kappa = ((0.01667*accFact-0.1644)*accFact+1.485)*accFact+0.3796;
        end

        sqrtAlpha = (1.0 + kappa * (1.0-sqrt(Tr)));
        alpha   = sqrtAlpha * sqrtAlpha ;
        aa      = 0.457235 *alpha* RR * RR * tCrit * tCrit / pCrit;
        bb      = 0.077796 * RR * tCrit / pCrit;

        RT = RR * TinK;
        A_ofP = aa*pres(i)/(RT*RT);
        B_ofP = bb * pres(i)/(RT);

        a2 = -(1.0d0 - B_ofP);
        a1 =  A_ofP - 3.0d0*B_ofP*B_ofP - 2.0d0*B_ofP;
        a0 = -B_ofP*( A_ofP - B_ofP*(1.0d0 + B_ofP) );

        [myRoots, numRoots] = rootsOfCubicEqn( a2, a1, a0 );

        if (numRoots==1)
            %Only one root
            z(i) = myRoots;     
        else
            %Selection of Root which gives lower Gibb's energy
            dG = RT*( (myRoots(2)-myRoots(1)) - log((myRoots(2)-B_ofP) ...
                 /(myRoots(1)-B_ofP))- A_ofP/(2.0*sqrt(2)*B_ofP) ...
                 *log((myRoots(2)+2.414*B_ofP)*(myRoots(1)-0.414*B_ofP) ...
                 /((myRoots(2)-0.414*B_ofP)*(myRoots(1)+2.414*B_ofP))) );
            if (dG > 0)
                z(i) = myRoots(1);
            else
                z(i) = myRoots(2);
            end
        end
    end
%     v = z*RT/pres;
%     
%     %compute fugacity of pure component
%     logPhi = z-1.0 -log(z-B_ofP) - A_ofP*log((z+2.414*B_ofP)/(z-0.414*B_ofP))...
%              /(2.0*B_ofP*sqrt(2));
end

