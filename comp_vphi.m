function vphi = comp_vphi(ts,te,p)

npl = length(p);      

vphi = zeros(1,npl*npl);

n = 1;

for i = 1:npl
    a = ts^(p(i));
    for j = 1:npl
        b = te^(p(j));
        vphi(n) = a*b;
        n = n+1;
    end
end
