function Vio = calVio(m,popg,poph,epsim)
Vio = zeros(m,1);
popg(popg<0) = 0;
if ~isempty(popg)
    Vio = Vio + sum(popg,2);
end
poph = abs(poph);
poph(poph<=epsim) = 0;
if ~isempty(poph)
    Vio = Vio + sum(poph,2);
end
end