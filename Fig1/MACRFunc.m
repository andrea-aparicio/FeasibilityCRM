function xdot = MACRFunc(x,params)
    A = params.A;
    r = params.r;
    xd = diag(x)*((A*x)+r);
%     xd= nan(1,length(A));
%     for i=1:length(xd)
%         xd(i)=x(i)*(1-(A(i,:)*x'))+D;
%     end

    xdot = xd;
end 
    

    