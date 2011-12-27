% Create surfaces from 3-D voxel data  

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%     Copyright (C) 2010  Ravi Soni www.hermesacademy.com ravi@hermesacademy.com
%
%     Program "Voxel Surface" is a free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     Program consists of the following files
%     run					Run a demo example
%     voxel_bnd_faces       Main function
%	  voxel_vtx             
%     renumbervertex
%     plotsurf              

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Input 
% voxels(nrow,ncol,nslice)  Ids of all voxels : nrow=ny ncol=nx nslice=nz
% sf(1,3)                   Voxel scale factors for x,y,z dimensions
%                           sf=[1 1 1] if problem does not have any dimensional scale
% offset(1,3)               Offset of crd of center of voxel(1,1,1) from global origin
%                           offset=[0 0 0] if global origin is the center of voxel(1,1,1)
% id                        Voxel id whose boundary is to be determined

% Output
% faces(nfaces,3)           Connectivity of boundary faces : vertices are numbered contigously from 1 to nvtx
% vertex(nvtx,dim)          Crd of boundary vertices : dim=3

% See fnc=voxel_bnd_faces for details

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

clear all
close hidden all

% Colors for plotting
red=[1.0 0.0 0.0];  
green=[0.5 1.0 0.5];  
cyan=[0.0 1.0 1.0];  
blue=[0.0 0.0 1.0]; 
yellow=[1.0 1.0 0.0];  
magenta=[1.0 0.0 1.0];
orange=[1.0 0.5 0.0];  
brown=[1.0 0.5 0.5]; 
purple=[0.5 0.5 1.0]; 
black=[0.0 0.0 0.0]; 
noclr='none';

% -------------------------------------------------------------------------
% Example 1 : CT Scan image with 5 voxel domains

load matf_imagesnifti  % image file

voxels=Images{1};  % voxels

voxid=1:5;  % voxel domain ids

% Scale factors for scaling each dimension
sf=[0.5508 0.5508 0.6000];
% 1 voxel along x axis=0.5508 units of length
% 1 voxel along y axis=0.5508 units of length
% 1 voxel along z axis=0.6 units of length
  
% Offset of crd of center of voxel(1,1,1) from global origin
offset=[0 0 0];  % center of voxel(1,1,1) = global origin

facecolor=[  % color of each domain
red
green
cyan
blue
yellow
]; 

% edgecolor='none';
edgecolor=[0 0 0];  % black

name='Surface';

fighndl=[];

% -------------------------------------------------------------------------
% Create surfaces for each voxel domain
% In this example faces and vertex are overwritten
% Store them in a cell array or structure if it is required to process them

for id=voxid
    id
    
    % Create surface from voxels
    [faces,vertex]=voxel_bnd_faces(voxels,sf,offset,id);
    
    % Plot surface
    [fighndl]=plotsurf(faces,vertex,facecolor(id,:),edgecolor,name,fighndl);  
        
end
