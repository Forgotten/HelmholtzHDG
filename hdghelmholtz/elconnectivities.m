function elcon = elconnectivities(nps,t2f)
% computed the element connectivities
%  nps: number of points in the edges of the elements
%  t2f triangle to faces matrix

% number of elements
ne  = size(t2f,1);

elcon = zeros(3*nps,ne);                % Element Connectivities
[il,jl] = find(t2f>0);

for i=1:length(il)
   f = t2f(il(i),jl(i));
   elcon((jl(i)-1)*nps+1:jl(i)*nps,il(i)) = (f-1)*nps+1:f*nps;
end

[il,jl] = find(t2f<0);

for i=1:length(il)
   f = -t2f(il(i),jl(i));
   elcon((jl(i)-1)*nps+1:jl(i)*nps,il(i)) = f*nps:-1:(f-1)*nps+1;
end

