    function fxi = fxi(f,x,i,xi)
    xx    = x;
    xx(i) = xi;
    fxi   = f(xx);
    end


