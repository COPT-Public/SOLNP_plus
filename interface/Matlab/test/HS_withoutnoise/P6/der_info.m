function der = der_info()
    der.gradf = @gradf;
    der.gradg = @gradg;
    der.hessf = @hessf;
    der.hessg = @hessg;
end

function gradf = gradf(x)
    gradf = [2*(x(1)-1);0];
end
function gradg = gradg(x)
    gradg = 10*[-2*x(1);1];
end
function hessf = hessf(x)
    hessf = [2,0;0,0];
end
function hessg = hessg(x)
    hessg = [-20,0;0,0];
end