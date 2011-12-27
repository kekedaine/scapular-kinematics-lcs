% Align axis of a base crd system to axis of a target crd system

% Axis1 of base crd system is aligned with Axis1 of target crd system, Axis2 with Axis 2 and Axis3 with Axis3.

% Rotations can be performed around either base crd system or target crd system : argument=reference

% Rotations can be performed in any axis order : argument=order
% Angles are output in the order of rotation.

% Note:
% If [thetax thetay thetaz] is a solution, 
% [thetax+pi thetay+pi thetaz+pi] & [thetax-pi thetay-pi thetaz-pi] are also solutions.

% Order may not yield a RH crd system, but a RH crd system is always enforced.
% Third axis of rotaxis and tgtaxis_ro are obtained from cross product of first two axis.
% If RH crd system is not enforced, results are incorrect.

% -------------------------------------------------------------------------
function [angles,angles_deg]=alignaxis(baseaxis,tgtaxis,reference,order)

% Input
% baseaxis(3,3)     Base crd system : rows are DC of each axis wrt global crd system
% tgtaxis(3,3)      Target crd system : rows are DC of each axis wrt global crd system
% reference         Crd system to be used as reference : 1=base  ~1=target
% order(3)          Order of axis for rotation : 1, 2 or 3

% Output
% angles(1,3)       Rotation angles around the 3 reference axis in radians
% angles_deg(1,3)   Rotation angles around the 3 reference axis in degrees

% Local
% rotaxis(3,3)      Rotated baseaxis : updated after each rotation
%                   Initially rotaxis=reordered baseaxis, finally rotaxis=tgtaxis_ro
% tgtaxis_ro(3,3)   Reordered target crd system

% -------------------------------------------------------------------------

% Define defaults
value=exist('reference','var');
if (value==0 | isempty(reference))
    % reference is undefined : use base crd system as default 
    reference=1;
end
    
value=exist('order','var');
if (value==0 | isempty(order))
    % order is undefined : use axis1,axis2,axis3 as default
    order=[1 2 3];
end

% Reorder axis
for i=1:2
    rotaxis(i,:)=baseaxis(order(i),:);
    tgtaxis_ro(i,:)=tgtaxis(order(i),:);
end
rotaxis(3,:)=cross(rotaxis(1,:),rotaxis(2,:));  % ensure RH crd system
tgtaxis_ro(3,:)=cross(tgtaxis_ro(1,:),tgtaxis_ro(2,:));  % ensure RH crd system
    
% Select reference crd system
if (reference==1)
    % Use base crd system as reference : rotations are about rotated baseaxis=rotaxis

    % ---------------------------------------------------------------------
    % Rotate about rotaxis1 (if required) such that rotaxis2 is in reordered target 1-2 plane
    
    % Find the vector in reordered target 1-2 plane that is normal to rotaxis1
    % This is the vector that rotaxis2 will be rotated to when rotated around rotaxis1
    vect12=cross(tgtaxis_ro(3,:),rotaxis(1,:));  % vector normal to tgtaxis_ro3 & rotaxis1
    if (norm(vect12)~=0)
        % Vectors are not parallel: rotate about rotaxis1 such that rotaxis2 is in reordered target 1-2 plane
        vect12=vect12/norm(vect12);  % unit vector normal to tgtaxis_ro3 & rotaxis1
        vec=rotaxis;
        vec1=rotaxis(2,:);
        vec2=vect12;
        pivot=[0 0 0];
            
        [rotaxis,normal1,theta1,theta1_deg]=rotvec(vec,vec1,vec2,pivot);
        
        dotprod=dot(rotaxis(1,:),normal1);
        if (dotprod<0)
            % Account for actual rotation about axis opposite to rotaxis1
            theta1=-theta1;
            theta1_deg=-theta1_deg;
        end
    else
        % Vectors are parallel : no rotation required about rotaxis1 as baseaxis2,baseaxis3 are in reordered target 1-2 plane
        rotaxis=rotaxis;
        theta1=0;
        theta1_deg=0;
    end   
    
    % ---------------------------------------------------------------------
    % Rotate about rotaxis2 such that rotaxis3=tgtaxis_ro3 (and rotaxis1 is in reordered target 1-2 plane) 
    vec=rotaxis;
    vec1=rotaxis(3,:);
    vec2=tgtaxis_ro(3,:);
    pivot=[0 0 0];
        
    [rotaxis,normal2,theta2,theta2_deg]=rotvec(vec,vec1,vec2,pivot);
    
    dotprod=dot(rotaxis(2,:),normal2);
    if (dotprod<0)
        % Account for actual rotation about axis opposite to rotaxis2
        theta2=-theta2;
        theta2_deg=-theta2_deg;
    end    

    % ---------------------------------------------------------------------
    % Rotate about rotaxis3 such that rotaxis1=tgtaxis_ro1 (and rotaxis2=tgtaxis_ro2)
    vec=rotaxis;
    vec1=rotaxis(1,:);
    vec2=tgtaxis_ro(1,:);
    pivot=[0 0 0];
      
    [rotaxis,normal3,theta3,theta3_deg]=rotvec(vec,vec1,vec2,pivot);
    
    dotprod=dot(rotaxis(3,:),normal3);
    if (dotprod<0)
        % Account for actual rotation about axis opposite to rotaxis3
        theta3=-theta3;
        theta3_deg=-theta3_deg;
    end     

else
    % Use target crd system as reference : rotations are about fixed tgtaxis 
    
    % ---------------------------------------------------------------------
    % Rotate about tgtaxis_ro1 (if required) such that rotaxis3 is in reordered target 1-3 plane
    
    % Project rotaxis3 onto reordered target 2-3 plane
    pnt=rotaxis(3,:);
    pntonplane=[0 0 0];
    normal=tgtaxis_ro(1,:);
    [dist,vec23,vecproj,location]=pntplanedist(pnt,pntonplane,normal);
    
    if (norm(vec23)~=0)
        % Rotate about tgtaxis1 such that rotaxis3 is in reordered target 1-3 plane
        % Accomplished by rotating vec23 onto tgtaxis_ro
        vec=rotaxis;
        vec1=vec23;
        vec2=tgtaxis_ro(3,:);
        pivot=[0 0 0];
        
        [rotaxis,normal1,theta1,theta1_deg]=rotvec(vec,vec1,vec2,pivot); 
        
        dotprod=dot(rotaxis(1,:),normal1);
        if (dotprod<0)
            % Account for actual rotation about axis opposite to rotaxis1
            theta1=-theta1;
            theta1_deg=-theta1_deg;
        end
    else
        % rotaxis3 is along tgtaxis_ro1 ; no rotation required around tgtaxis_ro1 
        rotaxis=rotaxis;
        theta1=0;
        theta1_deg=0;        
    end
    
    % ---------------------------------------------------------------------
    % Rotate about tgtaxis_ro2 such that rotaxis3=tgtaxis_ro3 (and rotaxis1-2 plane=tgtaxis_ro1-2 plane)
    vec=rotaxis;
    vec1=rotaxis(3,:);
    vec2=tgtaxis_ro(3,:);
    pivot=[0 0 0];
    
    [rotaxis,normal2,theta2,theta2_deg]=rotvec(vec,vec1,vec2,pivot); 

    dotprod=dot(rotaxis(2,:),normal2);
    if (dotprod<0)
        % Account for actual rotation about axis opposite to rotaxis1
        theta2=-theta2;
        theta2_deg=-theta2_deg;
    end
    
    % ---------------------------------------------------------------------
    % Rotate about tgtaxis_ro3 such that rotaxis1=tgtaxis_ro1 (and rotaxis2=tgtaxis_ro2)
    vec=rotaxis;
    vec1=rotaxis(1,:);
    vec2=tgtaxis_ro(1,:);
    pivot=[0 0 0];
    
    [rotaxis,normal3,theta3,theta3_deg]=rotvec(vec,vec1,vec2,pivot); 

    dotprod=dot(rotaxis(3,:),normal3);
    if (dotprod<0)
        % Account for actual rotation about axis opposite to rotaxis1
        theta3=-theta3;
        theta3_deg=-theta3_deg;
    end
    
end  % if (reference==1)
rotaxis;
tgtaxis_ro;

% Verify
[order_sort,indx]=sort(order,'ascend');
baseaxis_aligned=rotaxis(indx,:)
tgtaxis
tol=1e-6;
prod1=cross(baseaxis_aligned(1,:),tgtaxis(1,:));
prod2=cross(baseaxis_aligned(2,:),tgtaxis(2,:));
prod3=cross(baseaxis_aligned(3,:),tgtaxis(3,:));  % axis3 are anti-parallel if order does not yield RH crd system
if (norm(prod1)>tol | norm(prod2)>tol | norm(prod3)>tol)
    disp('Axis not aligned within tolerance')
end

% Output arguments
angles=[theta1,theta2,theta3];
angles_deg=[theta1_deg,theta2_deg,theta3_deg];

% -------------------------------------------------------------------------
% End of Function