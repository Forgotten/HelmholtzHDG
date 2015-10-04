function [f,t2f]=mkt2f(t)
%MKT2F Compute Face Connectivity and Triangle to Face Connetivity.
%   [F,T2F]=MKT2F(T)
%
%      T:         Triangle indices (NT,3)
%      F:         Face connectivity (NF,4) (for boundary edges F(:,4)=0)
%      T2F:       Triangle to Face Connectivity (NT,3)
%
%   See also: MKT2T.
t2t = mkt2t(t);
nb = sum(sum(t2t <= 0));
f = zeros((3*size(t,1)+nb)/2,4);
t2f = zeros(size(t));
jf = 0;
for i=1:size(t,1)
    for j = 1:3
        if t2t(i,j) > i || t2t(i,j) <=0
            ie = t2t(i,j);
            jf = jf + 1;
            f(jf,1) = t(i,mod(j,3)+1);
            f(jf,2) = t(i,mod(j+1,3)+1);
            f(jf,3) = i;
            f(jf,4) = ie;
            t2f(i,j) = jf;
            if ie > 0
               if (f(jf,1) > f(jf,2)) % Orient interior faces mesh.f(:,1) < mesh.f(:,2)
                   f(jf,:) = f(jf,[2,1,4,3]);
                   t2f(i,j) = -jf;
                   t2f(ie,find(t(ie,:)==(sum(t(ie,:))-sum(f(jf,1:2))))) = jf;
               else
                   t2f(ie,find(t(ie,:)==(sum(t(ie,:))-sum(f(jf,1:2))))) = -jf;
               end
            end
        end
    end
end
clear t2t;

% Reorder faces - First interior then boundary

nf = size(f,1);
ne = size(t2f,1);
[a,mp] = sort(f(:,4) == 0);
f = f(mp,:);

a = zeros(nf,1);
a(mp) = (1:nf)';
t2f = reshape(t2f,3*ne,1);
ii = find(t2f > 0);
t2f(ii) = a(t2f(ii));
ii = find(t2f < 0);
t2f(ii) = -a(-t2f(ii));
t2f = reshape(t2f,ne,3);





