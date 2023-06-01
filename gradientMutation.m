function [y]=gradientMutation(x,gg,hh,funcno,func)
    % preset parameters of forward differences approach for Jacobian approximation
    par = Cal_par(funcno);
    eta     = 1e-4;                            
    n       = length(x);
    Dd      = eye(n);
    dx      = repmat(x,1,n)+eta.*Dd;
    for i=1:n
        dx(dx(:,i)'>par.xmax,i) = par.xmax(dx(:,i)'>par.xmax);
        [fff, gv, hv] = func(dx(:,i)',funcno);
	    dg= [gv';hv'];
	    dCx(:,i)  = dg;
    end
    
    %create vector of coanstraint-wise constraint violations
    deltaG  = max(0,gg);
    Cx= [gg;hh];
    
    % approaximate Jacobian
    nabC    = 1/eta.*( dCx - repmat(Cx,1,n));
    delC    = [deltaG;hh];
    
    % compute Moore-Penrose inverse of nabC
    try
        inv_nabC= pinv(nabC,1e-12);  
        deltaX  = -inv_nabC*delC;
    catch
        deltaX  = 0;
    end
    
    % repair the infeasible candidate solution x
    y       = (x+deltaX);
end
