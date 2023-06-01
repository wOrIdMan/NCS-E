function RelaxedConV = RelaxConV(m,g,h,epsim,v)
RelaxedConV = zeros(m,1);
g(g<v) = 0;
if ~isempty(g)
    RelaxedConV = RelaxedConV + sum(g,2);
end
h = abs(h);
h(h<=epsim+v) = 0;
if ~isempty(h)
    RelaxedConV = RelaxedConV + sum(h,2);
end
end