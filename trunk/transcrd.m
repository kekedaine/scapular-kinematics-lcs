% Function for transforming coordinates from base crd system to target crd system

% Only rotation is accounted for in this function, so translate crd externally as required
% This has been done so that the function can be used for reverse transformation as well

% To transform crd from target to base crd system use transpose of the matrix axis

% -------------------------------------------------------------------------

function [crd]=transcrd(crdbase,axis)

% Input
% crdbase (n,3)            Crds w.r.t base crd system
% axis(3,3)                Transformation matrix from base to target : col are DC of axis of target crd system w.r.t base crd system
%                             x-axis y-axis z-axis
%                          l  lx     ly     lz
%                          m  mx     my     mz
%                          n  nx     ny     nz

% Output
% crd(n,3)                 Crd w.r.t target crd system

% -------------------------------------------------------------------------

% (x,y,z)=base crd of point

% target crd i= projection of point vector on target axis i = dot(point crd,dc of target axis i)
% target crd x= x*lx+y*mx+z*nx  
% target crd y= x*ly+y*my+z*ny 
% target crd z= x*lz+y*mz+z*nz 

[n,d]=size(crdbase);

crd=crdbase*axis;   % project position vector onto target axis

% -------------------------------------------------------------------------
% End of Function