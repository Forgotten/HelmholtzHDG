function u = initu(mesh,app,param)
%INITU Initialize vector of unknowns
%    U=INITU(MESH,APP,PARAM)
%
%    mesh:             Mesh structure
%    app:              Application structure
%    param{APP.NC}:    Cell array containing
%                      When param{i} is a double U(:,APP.NC,:) = param{i}
%                      When param{i} is a pointer to a function,
%                                    U(:,APP.NC,:) = param{i}(mesh.dgnodes)
%    U(NPL,APP.NC,NT): Scalar function to be plotted
%
u = zeros(size(mesh.dgnodes,1),app.nc,size(mesh.dgnodes,3));
for i=1:app.nc
    if isa(param{i},'double')
        u(:,i,:) = param{i};
    else
        u(:,i,:) = param{i}(mesh.dgnodes);
    end
end

end