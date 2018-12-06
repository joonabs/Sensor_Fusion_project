% Calculates the "optimal" step length using golden sectioning
function minimum = golden_sectioning(f, R, L)
    a1 = L + (1 - (-1 + sqrt(5))/2)*abs(R-L);
    a2 = L + ((-1 + sqrt(5))/2)*abs(R-L);

    while abs(R - L) > 0.01
      if f(a1) < f(a2)
          R = a2;
          a2 = a1;
          a1 = L + (1 - (-1 + sqrt(5))/2)*abs(R-L);
      else
          L = a1;
          a1 = a2;
          a2 = L + ((-1 + sqrt(5))/2)*abs(R-L);
      end
    end
    
    if f(a1) < f(a2)
        minimum = a1;
    else
        minimum = a2;
    end
end