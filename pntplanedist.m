% Normal distances of a point to planes with a common point and given normals
% Projected points, projected vectors and location of point wrt planes are also computed

% -------------------------------------------------------------------------

function [dist,pntproj,vecproj,location]=pntplanedist(pnt,pntonplanes,normals)

% Input
% pnt(1,dim)                Coordinate of point
% pntonplanes(1,dim)        Coordinate of common point on planes
% normals(nfc,dim)          Normal vectors to each plane

% Output
% dist(nfc,1)               Normal distance of point to each plane
% pntproj(nfc,3)            Projection of point on each plane
% vecproj(nfc,3)            Unit vector from pntproj to pnt
% location(nfc,1)           Location flag for point
%                           0=on plane 
%                           1=towards normal
%                           -1=opposite to normal
% Local
% tol                       Tolerance for determining location=0
%                           Use very small tolerance (<<min edge length) to include only points exactly on the plane
%                           Use small tolerance (~min edge length) to include vertices of faces that intersect the plane

% -------------------------------------------------------------------------

[nfc,dim]=size(normals);

tol=0.0001; 
% tol=0.7;  

for i=1:nfc
    % Compute normal distance to each plane
    
    % Plane equation general form : ax+by+cz+d=0  [a b c]=normal  [x y z]=point onplane
    normali=normals(i,:);
    pntonplanes;
    d=-dot(normali,pntonplanes);
    
    % Projection of point on plane
    s=square(norm(normali));  % normals are normalized so s should be 1
    k=-(dot(normali,pnt)+d)/s;
    pntproj(i,:)=pnt+normali*k;
    
    % Normal distance
    del=pnt-pntproj(i,:);
    dist(i,1)=norm(del);
    
    % Unit projected vector
    if (dist(i,1)<tol)
        % Point on plane
        vecproj(i,:)=[0 0 0];    
        location(i,1)=0;
    else
        % Exterior or interior point
        vecproj(i,:)=del/dist(i,1); 
        dotprod=dot(vecproj(i,:),normali);
        if (dotprod>0)
            location(i,1)=1;  % exterior point
        else
            location(i,1)=-1; % interior point
        end
    end
end
pntproj;
location;

% -------------------------------------------------------------------------
% End of Function