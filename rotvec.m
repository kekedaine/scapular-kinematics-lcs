% Function for rotating vector using two reference vectors

% Axis of rotation  = normal to plane of reference vectors and passes thru point=pivot
% Angle of rotation = angle between reference vectors (theta)

% A target crd system is constructed such that 
% xt=ref vec1  
% zt=common normal of reference vectors
% yt=zt x xt

% A rotated crd system is constructed such that 
% zr=zt
% xr=xt rotated by -theta around zr=zt
% yr=yt rotated by -theta around zr=zt

% Rotation is achieved in 5 steps:
% 1. Translate vector by delta=-pivot
% 2. Transform vector to target crd system
% 3. Rotate vector in target crd system (thru crd transformation - see below)
% 4. Transform rotated vector to base crd system
% 5. Translate rotatedvector by delta=+pivot

% If vector=(r,phi,z) in target crd system, then rotated vector=(r,phi+theta,z)
% Consider a crd system that is obtained by rotating target crd system by -theta around z
% Then crd of unrotated vector about the rotated crd system=(r,phi+theta,z)
% Thus we can achieve rotation by theta by transforming to a crd system rotated by -theta
%       vectr = Ttr(-theta)*vect
%             = Trt(theta)*vect   reverse transformation is inverse operation
%             = axisr*vect          
%
%       vectr = crd of rotated vector in target crd system
%       vect  = crd of unrotated vector in target crd system
%       Ttr   = transformation matrix from target to rotated crd system 
%       Trt   = transformation matrix from rotated to target crd system 
%             = axisr
 
% Note:
% If vector represents a point and the point has to be rotated and moved to another point,
% then first use this function to rotate and then translate the point to the new position extrernally

% -------------------------------------------------------------------------

function [vecr,normal,theta,theta_deg]=rotvec(vec,vec1,vec2,pivot)

% Input
% vec(n,3)          Vectors to be rotated  
% vec1(1,3)         Reference vector 1 
% vec2(1,3)         Reference vector 2 
% delta(1,3)        Offset of axis of rotation wrt tail of vec
% pivot(1,3)        Pivot point for rotation

% Output
% vecr(n,3)         Rotated vector
% normal(1,3)       Common normal of reference vectors
% theta             Angle between reference vectors in radians
% theta_deg         Angle between reference vectors in degrees

% Local
% axistb(3,3)       Transformation matrix from base to target : dc of axis of target crd system w.r.t base crd system 
%                   Target crd system : xt=vec1  zt=common normal  yt=zt x xt
%                       l  m  n
%                   xt  lx mx nx
%                   yt  ly my ny
%                   zt  lz mz nz

% axisr(3,3)        Rotation matrix about zt = transformation matrix from rotated to target crd system
%                   l                m                 n
%                   cos(angle)       sin(angle)        0  
%                   cos(angle+pi/2)  sin(angle+pi/2)   0
%                   0                0                 1  

% axisbt(3,3)       Transformation matrix from target to base = Transpose of axistb 

% -------------------------------------------------------------------------

[n,d]=size(vec);

normal=cross(vec1,vec2);  % common normal
dotprod12=dot(vec1,vec2);

norm_normal=norm(normal);

if (norm_normal~=0)
    % vec1 & vec2 are not collinear : rotate vec

    % Angle between vectors (do not use sin as it will always return acute angle)
    cos_theta=dotprod12/norm(vec1)/norm(vec2);
    theta=real(acos(cos_theta));  % small imaginary part can result if cos_theta mag is slightly greater than 1 due to roundoff errors
    
    theta_deg=(180/pi)*theta;
    normal=normal/norm_normal;
    
    xt=vec1/norm(vec1);
    zt=normal;
    yt=cross(zt,xt);

    axistb=[xt;yt;zt]; % Transformation matrix from target to base : col are DC of base axis wrt target axis (rows are DC of target axis wrt base axis)
    axisbt=axistb.';   % Transformation matrix from base to target : col are DC of target axis wrt base axis 

    % Transform vector to target crd system
    for i=1:n
        vec(i,:)=vec(i,:)-pivot;
    end
    vec;
    [vect]=transcrd(vec,axisbt); % transform to target crd system

    % Rotate vector by angle=theta about z-axis of target crd system (common normal) 
    % Achieved by transforming to a crd system rotated about z-axis of target crd system by -theta
    % This is equivalent to transforming from rotated crd system to target crd system
    % rotation(theta)=Ttr(-theta)=Trt(theta)=axisr
    angle=theta;   
    %     l                 m                 n
    axisr=[cos(angle)       sin(angle)        0;  % Trt=Transformation matrix from rotated to target crd system  
          cos(angle+pi/2)   sin(angle+pi/2)   0;  %     Col are DC of target axis wrt rotated axis 
          0                 0                 1]; %     Rows are DC of rotated axis wrt target axis    

    [vectr]=transcrd(vect,axisr);

    % Transform rotated vector to global crd system
    [vecr]=transcrd(vectr,axistb); % transform to base crd system
    for i=1:n
        vecr(i,:)=vecr(i,:)+pivot;
    end
    vecr;
else
    % vec1 & vec2 are collinear : do not rotate vec
    vecr=vec;
    theta=0;
    theta_deg=0;
end

% -------------------------------------------------------------------------
% End of Function