clear all
clc
con_output = [];
MAHD_output = [];
sca_ver_output = [];
pre_path = 'dat\t2\';
allfiles = dir([pre_path 's15c1.nii']);

% BEGIN PARAMETER SETUP
%% Global
run_conoid_analysis = 1; % turn on(1)/off(0) conoid_analysis
run_MAHD_analysis = 1; % turn on(1)/off(0) MAHD_analysis
run_sca_ver_analysis = 1; % turn on(1)/off(0) sca_ver_analysis

first_cut_plane_by_ABC = 0; % = 0 to indicate not used
first_cut_plane_by_PCA = 0; % = 0 to indicate not used
first_cut_plane_by_AD = 1; % = 0 to indicate not used
%% Conoid parameters
max_conoid_sv_precision = 1; % subvoxel precsion
PCA_axis = 3; % which PCA axis is chosen to idetify the cutting angle (valid values are 1,2,3)
acr_trim_rate = 0.333; % rate of to-be-kept volume when trimming acromion (the orginal document specified this param as 0.33)
render_con_output = 1; % = 0 to turn off rendering (including .fig)

%% MAHD parameters
max_MAHD_sv_precision = 1; % subvoxel precsion

%% Scapula - Vertebra parameters
invert_vertebra_y_axis = 0; % = 1 to invert the direction of vector y, = 0 to use default direction of vector y

% set one of the following flag to 1 and the rest to 0 to select the option
use_ABC_to_define_Scapula_Vertebra_LCS = 0;
use_DEF_to_define_Scapula_Vertebra_LCS = 0;
use_DEF_TPI_TPC_SP_to_define_Scapula_Vertebra_LCS = 1;
% END PARAMETER SETUP

%% Common processing
for file_id = 1:size(allfiles,1),
close all;
filename = [pre_path allfiles(file_id).name],

total_time = tic;

acr_trim_rate = max(min(1,acr_trim_rate),0.333);
PCA_axis = mod(PCA_axis-1,3)+1;

PRECISION = sqrt(3);

AC_ID = 1; % acromion label = 1
HU_ID = 2; % humerus label = 2
SC_ID = 3; % scapula label = 3
VE_ID = 4; % vertebra label = 4
TA_ID = 6; % trimmed acromion label = 5
SA_ID = 11; % label of subarcomial space (the part of the conoid outside the acromion and humerus) = 11
CA_ID = TA_ID+SA_ID; % label of the part of the conoid inside the trimmed acromion = 17
CH_ID = HU_ID+SA_ID; % label of the part of the conoid inside the humerus = 13
A_ID = 10; % label of point A = 10
B_ID = 7; % label of point B = 7
C_ID = 9; % label of point C = 9
D_ID = 7; % label of point D = 7
E_ID = 9; % label of point E = 9
F_ID = 12; % label of point F = 12
TPI_ID = 14; % transverse process ipsilateral
TPC_ID = 15; % transverse process contralateral
SP_ID = 16; % spinous process

if (use_DEF_to_define_Scapula_Vertebra_LCS&&(B_ID==D_ID||C_ID==E_ID))
    first_cut_plane_by_ABC = 0; % if no B or C, we cannot use ABC to trim the acromion
    first_cut_plane_by_AD = 1;
    first_cut_plane_by_PCA = 0;
end

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

facecolor = [red; green; cyan; blue; yellow; magenta; orange; brown; purple; black];

scan = load_nii(filename);
os = size(scan.img);
scale_factor = scan.hdr.dime.pixdim(2:4);
cvt_precision = 1/min(scale_factor); % so no voxel disappears after space transform
s = round(os .* scale_factor * cvt_precision);
n_voxel = s(1) * s(2) * s(3);
imgo = uint8(scan.img);
clear scan;
imgd = zeros(s(1),s(2),s(3),'uint8');

% convert metric space
lut1 = double(1:s(1));
%lut1 = ceil(lut1 * double(os(1)) / double(s(1)));
lut1 = round((lut1-1) * (double(os(1)) / double(s(1))) + 0.5); 
lut2 = double(1:s(2));
%lut2 = ceil(lut2 * double(os(2)) / double(s(2)));
lut2 = round((lut2-1) * (double(os(2)) / double(s(2))) + 0.5);
lut3 = double(1:s(3));
%lut3 = ceil(lut3 * double(os(3)) / double(s(3)));
lut3 = round((lut3-1) * (double(os(3)) / double(s(3))) + 0.5);

imgd(:,:,:) = imgo(lut1,lut2,lut3);
clear lut1 lut2 lut3 imgo;

step = 1; % speed-precision control
% step 0: identify surface voxels using neighbour-6
tic
fprintf('Indexing surface voxels...');
surf_voxels = zeros(n_voxel/27,4);
hum_voxels = zeros(n_voxel/27,3);
acr_voxels = zeros(n_voxel/27,3);

id_s = 0;
id_h = 0;
id_a = 0;
for x=1:step:s(1);
    for y=1:step:s(2);
        for z=1:step:s(3);
            v = imgd(x,y,z);
            if (v ~= 0)
                if (v == 2)
                    id_h = id_h + 1;
                    hum_voxels(id_h,:) = [x y z];
                else if (v == 1)
                    id_a = id_a + 1;
                    acr_voxels(id_a,:) = [x y z];
                    end
                end
                if (x==1 || x==s(1) || y==1 || y==s(2) || z==1 || z==s(3) || ((v ~= imgd(x,y,z-1)) || (v ~= imgd(x,y,z+1)) || (v ~= imgd(x,y-1,z)) || (v ~= imgd(x,y+1,z)) || (v ~= imgd(x-1,y,z)) || (v ~= imgd(x+1,y,z))))
                    id_s = id_s + 1;
                    surf_voxels(id_s,:) = [x y z double(v)];
                end
            end
        end
    end
end

surf_voxels = surf_voxels(1:id_s,:);
hum_voxels = hum_voxels(1:id_h,:);
acr_voxels = acr_voxels(1:id_a,:);
% conoid apex
apex = (sum(hum_voxels)/size(hum_voxels,1));
fprintf('done in %4.4f s\n',toc);

% step 1 trim the size of the segmented acromion
% angle of cutting plane
tic
fprintf('Trimming the acromion...');
pos_acr_sur_id = find(surf_voxels(:,4) == A_ID);
if (size(pos_acr_sur_id, 1) < 1)
    error('\n Most posterior acromial point not found!\n')
end
pos_acr_sur_id = pos_acr_sur_id(1);

sca_sur_id = find(surf_voxels(:,4) == SC_ID);
acr_sur_id = [find(surf_voxels(:,4) == AC_ID); pos_acr_sur_id];
pos_acr_sur_id = size(acr_sur_id,1);

sca_acr_sur_id = [acr_sur_id; sca_sur_id];
[SAPCA,sc] = princomp(surf_voxels(sca_acr_sur_id,1:3));

pos_SAPCA_Z = sc(pos_acr_sur_id,3);

% use PCA to find most anterior point
% [temp,ant_acr_sur_id] = max(abs(sc(1:pos_acr_sur_id,3) - pos_SAPCA_Z));
% if (size(ant_acr_sur_id,1) < 1)
%     error('\n Most anterior acromial point not found!\n')
% end
% ant_acr_sur_id = ant_acr_sur_id(1);
% 
% % debug ant and pos point
% if (render_con_output == 1)
% [x y z] = bresenham_line3d(surf_voxels(acr_sur_id(ant_acr_sur_id),1:3), surf_voxels(acr_sur_id(pos_acr_sur_id),1:3),2);
% figure(10)
% hold on
% plot3(x,y,z);
% end

%% Scapula - Vertebra analysis
if (run_sca_ver_analysis == 1)
    pA = surf_voxels(surf_voxels(:,4)==A_ID,1:3); pA = pA(1,:);
    
    sca_sur = surf_voxels(sca_sur_id,1:3);
    sca_or = pA;
    
    if (~use_DEF_TPI_TPC_SP_to_define_Scapula_Vertebra_LCS)
    ver_sur = surf_voxels(surf_voxels(:,4)==VE_ID,1:3);
    ver_or = mean(ver_sur);
    VPCA = princomp(ver_sur);

    %rearange VPCA so that xyz = [secondary tetiary primary]
    VLCS = [VPCA(:,2) VPCA(:,3) VPCA(:,1)];
    end
        
    if (use_ABC_to_define_Scapula_Vertebra_LCS)

        pB = surf_voxels(surf_voxels(:,4)==B_ID,1:3); pB = pB(1,:);
        pC = surf_voxels(surf_voxels(:,4)==C_ID,1:3); pC = pC(1,:);
    
        if (invert_vertebra_y_axis == 1), input_sign_ver_y = -1; else input_sign_ver_y = 1; end

        pC_sign = (pC - ver_or) * VLCS(:,1); if (pC_sign > 0);  VLCS(:,1) = - VLCS(:,1); end% #1 axis points positive away from C
        pC_sign = (pC - ver_or) * VLCS(:,2); if (pC_sign > 0);  VLCS(:,2) = - VLCS(:,2); end
        VLCS(:,2) = VLCS(:,2) * input_sign_ver_y;% #2 axis points positive away from C, we also consider the input sign (if any)
        pC_sign = (pC - ver_or) * VLCS(:,3); if (pC_sign < 0);  VLCS(:,3) = - VLCS(:,3); end% #3 axis points C positive

        figure(10);
        arrow3d(ver_or,ver_or+VLCS(:,1)'*64,red,red);
        arrow3d(ver_or,ver_or+VLCS(:,2)'*64,green,green);
        arrow3d(ver_or,ver_or+VLCS(:,3)'*64,blue,blue);

        sca_az = (pA - pB)'; sca_az = sca_az / norm(sca_az);
        sca_ax = null([sca_az'; (pA - pC)]); % perpendicular to ABC
        ver_sign = (ver_or - sca_or) * sca_ax; if (ver_sign < 0); sca_ax = -sca_ax; end% point anteriorly
        sca_ay = null([sca_az'; sca_ax']); if (sca_ay*(pB - pC) < 0); sca_ay = -sca_ay; end% perpendicular to az and ax and point superiorly
        SLCS = [sca_ax sca_ay sca_az];
    
    end
    
    if (use_DEF_to_define_Scapula_Vertebra_LCS)
        
        pD = surf_voxels(surf_voxels(:,4)==D_ID,1:3); pD = pD(1,:);
        pE = surf_voxels(surf_voxels(:,4)==E_ID,1:3); pE = pE(1,:);
        pF = surf_voxels(surf_voxels(:,4)==F_ID,1:3); pF = pF(1,:);
        
        if (invert_vertebra_y_axis == 1), input_sign_ver_y = -1; else input_sign_ver_y = 1; end

        % base VLCS on point F
        pF_sign = (pF - ver_or) * VLCS(:,1); if (pF_sign < 0);  VLCS(:,1) = - VLCS(:,1); end% #1 axis points F positive
        pF_sign = (pF - ver_or) * VLCS(:,2); if (pF_sign > 0);  VLCS(:,2) = - VLCS(:,2); end
        VLCS(:,2) = VLCS(:,2) * input_sign_ver_y;% #2 axis points positive away from F, we also consider the input sign (if any)
        pF_sign = (pF - ver_or) * VLCS(:,3); if (pF_sign < 0);  VLCS(:,3) = - VLCS(:,3); end% #3 axis points F positive

        figure(10);
        arrow3d(ver_or,ver_or+VLCS(:,1)'*64,red,red);
        arrow3d(ver_or,ver_or+VLCS(:,2)'*64,green,green);
        arrow3d(ver_or,ver_or+VLCS(:,3)'*64,blue,blue);

        % base SLCS on points A, D, E
        sca_ax = (pD - pA)'; sca_ax = sca_ax / norm(sca_ax); %  the x-axis is formed by line AD (positive points forward)
        sca_ay = null([sca_ax'; (pE - pA)]); % perpendicular to ADE
        if (sca_ay(3)<0); sca_ay = sca_ay*-1; end % y axis point superior positive, here we assume the GCS is xxS (z axis points superiorly)
        sca_az = null([sca_ax'; sca_ay']);
        ver_sign = (ver_or - sca_or) * sca_az; if (ver_sign > 0); sca_az = -sca_az; end% the z-axis is the common line perpendicular to both the y and x axes, with positive pointing laterally or away from vertebra
        SLCS = [sca_ax sca_ay sca_az];
        
    end
    
    if (use_DEF_TPI_TPC_SP_to_define_Scapula_Vertebra_LCS)
        
        pD = surf_voxels(surf_voxels(:,4)==D_ID,1:3); pD = pD(1,:);
        pE = surf_voxels(surf_voxels(:,4)==E_ID,1:3); pE = pE(1,:);
        pF = surf_voxels(surf_voxels(:,4)==F_ID,1:3); pF = pF(1,:);
        
        pTPI = surf_voxels(surf_voxels(:,4)==TPI_ID,1:3); pTPI = pTPI(1,:);
        pTPC = surf_voxels(surf_voxels(:,4)==TPC_ID,1:3); pTPC = pTPC(1,:);
        pSP = surf_voxels(surf_voxels(:,4)==SP_ID,1:3); pSP = pSP(1,:);
        
        if (invert_vertebra_y_axis == 1), input_sign_ver_y = -1; else input_sign_ver_y = 1; end

        % base VLCS on point TPI(14) TPC(15) SP(16)
        % Origin is the midpoint of of the line connecting points 14 and 15
        ver_or = (pTPI + pTPC)/2;
        % z-axis is the line connecting points 14 and 15
        VLCS(:,3) = (ver_or - pTPI)';
        VLCS(:,3) = VLCS(:,3) / norm(VLCS(:,3));
        pF_sign = (pF - ver_or) * VLCS(:,3); if (pF_sign < 0);  VLCS(:,3) = - VLCS(:,3); end% #3 axis points F positive
        % y-axis is the line perpendicular to the plane formed by connecting points 14, 15 and 16
        VLCS(:,2) = null([(pTPI-pTPC); (pTPC-pSP)]);
        pF_sign = (pF - ver_or) * VLCS(:,2); if (pF_sign < 0);  VLCS(:,2) = - VLCS(:,2); end% #2 axis points F positive
        % x-axis is the line perpendicular to both z-axis and y-axis
        VLCS(:,1) = null([VLCS(:,2)'; VLCS(:,3)']);
        pF_sign = (pF - ver_or) * VLCS(:,1); if (pF_sign < 0);  VLCS(:,1) = - VLCS(:,1); end% #1 axis points F positive

        figure(10);
        arrow3d(ver_or,ver_or+VLCS(:,1)'*64,red,red);
        arrow3d(ver_or,ver_or+VLCS(:,2)'*64,green,green);
        arrow3d(ver_or,ver_or+VLCS(:,3)'*64,blue,blue);

        % base SLCS on points A, D, E
        sca_ax = (pD - pA)'; sca_ax = sca_ax / norm(sca_ax); %  the x-axis is formed by line AD (positive points forward)
        sca_ay = null([sca_ax'; (pE - pA)]); % perpendicular to ADE
        if (sca_ay(3)<0); sca_ay = sca_ay*-1; end % y axis point superior positive, here we assume the GCS is xxS (z axis points superiorly)
        sca_az = null([sca_ax'; sca_ay']);
        ver_sign = (ver_or - sca_or) * sca_az; if (ver_sign > 0); sca_az = -sca_az; end% the z-axis is the common line perpendicular to both the y and x axes, with positive pointing laterally or away from vertebra
        SLCS = [sca_ax sca_ay sca_az];
        % special case: cutting plane is parallel to AE
        cut_norm = null([sca_ay'; (pE - pA)]);
        
        % debug D,E
        pD_cube = [pD; pD+[0 1 0]; pD+[1 0 0]; pD+[0 -1 0]; pD+[-1 0 0]];
        pE_cube = [pE; pE+[0 1 0]; pE+[1 0 0]; pE+[0 -1 0]; pE+[-1 0 0]];
        figure(10)
        hold on
        [t]=delaunay(pD_cube(:,1),pD_cube(:,2));
        trisurf(t,pD_cube(:,1),pD_cube(:,2),pD_cube(:,3),'facecolor',[1 0 0],'edgecolor',[0.8 0.8 0.8]); % pD is RED
        arrow3d(pA,pD,red,red);
        hold on
        [t]=delaunay(pE_cube(:,1),pE_cube(:,2));
        trisurf(t,pE_cube(:,1),pE_cube(:,2),pE_cube(:,3),'facecolor',[0 1 0],'edgecolor',[0.8 0.8 0.8]); % pE is GREEN
        arrow3d(pA,pE,green,green);
    end

    figure(10);
    arrow3d(sca_or,sca_or+sca_ax'*64,red,red);
    arrow3d(sca_or,sca_or+sca_ay'*64,green,green);
    arrow3d(sca_or,sca_or+sca_az'*64,blue,blue);

    % find the rotation matrix to align SLCS axes to VLCS axes
    mobile_dcm = SLCS';
    static_dcm = VLCS';
    % find the rotation matrix to align SLCS axes to VLCS axes
    R = static_dcm*mobile_dcm';
    [Yr Xr Zr]=dcm2angle(R,'YXZ');
    S2V_angles_deg = [Yr Xr Zr]/pi*180;

    figure(9)
    axis equal
    arrow3d([0 0 0],VLCS(:,1)',red,black);
    arrow3d([0 0 0],VLCS(:,2)',green,black);
    arrow3d([0 0 0],VLCS(:,3)',blue,black);
    %draw xy plane
    X = [0 VLCS(1,[1 2])];
    Y = [0 VLCS(2,[1 2])];
    Z = [0 VLCS(3,[1 2])];
    C = [1 1 0];
    patch(X,Y,Z,C,'facealpha',0.25);

    arrow3d([0 0 0],SLCS(:,1)',red,red);
    arrow3d([0 0 0],SLCS(:,2)',green,green);
    arrow3d([0 0 0],SLCS(:,3)',blue,blue);
    %draw XZ plane
    X = [0 SLCS(1,[1 3])];
    Y = [0 SLCS(2,[1 3])];
    Z = [0 SLCS(3,[1 3])];
    C = [0 1 1];
    patch(X,Y,Z,C);

    % find the line where the 2 planes meet
    direction_vect = null([sca_ay'; VLCS(:,3)']);
    arrow3d([0 0 0],direction_vect',black,black);
    arrow3d([0 0 0],-direction_vect',black,black);
    if (first_cut_plane_by_ABC == 1)
        sca_az = (pA - pB)'; sca_az = sca_az / norm(sca_az);
        sca_ax = null([sca_az'; (pA - pC)]); % perpendicular to ABC
        ver_sign = (ver_or - sca_or) * sca_ax; if (ver_sign < 0); sca_ax = -sca_ax; end% point anteriorly
        sca_ay = null([sca_az'; sca_ax']); if (sca_ay*(pB - pC) < 0); sca_ay = -sca_ay; end% perpendicular to az and ax and point superiorly
        SAPCA = [sca_az sca_ay sca_ax]; % only use the third vector (sca_ax) so just need to make sure it is a unit vector perpendicular to first cutting plane
    else if (first_cut_plane_by_AD == 1)
            SAPCA = [cut_norm cut_norm cut_norm];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use AE to find most anterior point
[temp,ant_acr_sur_id] = max(surf_voxels(sca_acr_sur_id,1:3)*sca_az);
if (size(ant_acr_sur_id,1) < 1)
    error('\n Most anterior acromial point not found!\n')
end
ant_acr_sur_id = ant_acr_sur_id(1);

% debug ant and pos point
if (render_con_output == 1)
[x y z] = bresenham_line3d(surf_voxels(acr_sur_id(ant_acr_sur_id),1:3), surf_voxels(acr_sur_id(pos_acr_sur_id),1:3),2);
figure(10)
hold on
plot3(x,y,z);
end

% 1st cutting plane
cut_point = (surf_voxels(acr_sur_id(pos_acr_sur_id),1:3) * acr_trim_rate + (1-acr_trim_rate) * surf_voxels(acr_sur_id(ant_acr_sur_id),1:3)); % this part is used for PCA
cut_point_base = (surf_voxels(acr_sur_id(pos_acr_sur_id),1:3) * 0.334 + 0.666 * surf_voxels(acr_sur_id(ant_acr_sur_id),1:3)); % 1/3 from anterior

% debug cut_point1
if (render_con_output == 1)
cut_point_base_cube = [cut_point_base; cut_point_base+[0 1 0]; cut_point_base+[1 0 0]; cut_point_base+[0 -1 0]; cut_point_base+[-1 0 0]];
cut_point_cube = [cut_point; cut_point+[0 1 0]; cut_point+[1 0 0]; cut_point+[0 -1 0]; cut_point+[-1 0 0]];
figure(10)
hold on
[t]=delaunay(cut_point_cube(:,1),cut_point_cube(:,2));
trisurf(t,cut_point_cube(:,1),cut_point_cube(:,2),cut_point_cube(:,3),'facecolor',[1 0 0],'edgecolor',[0.8 0.8 0.8]);
hold on
[t]=delaunay(cut_point_base_cube(:,1),cut_point_base_cube(:,2));
trisurf(t,cut_point_base_cube(:,1),cut_point_base_cube(:,2),cut_point_base_cube(:,3),'facecolor',[0 1 0],'edgecolor',[0.8 0.8 0.8]);
end

d_pos = [surf_voxels(acr_sur_id(pos_acr_sur_id),1)-cut_point(1) surf_voxels(acr_sur_id(pos_acr_sur_id),2)-cut_point(2) surf_voxels(acr_sur_id(pos_acr_sur_id),3)-cut_point(3)] * SAPCA(:,3);
sign_pos = sign(d_pos);

% mark the surface of trimmed acromion
d = [surf_voxels(acr_sur_id,1)-cut_point_base(1) surf_voxels(acr_sur_id,2)-cut_point_base(2) surf_voxels(acr_sur_id,3)-cut_point_base(3)] * SAPCA(:,3);
trim_acr_sur_id = find(sign(d) ~= sign_pos);
surf_voxels(acr_sur_id(trim_acr_sur_id),4) = TA_ID;

% trim the acromion 1/3
d = [acr_voxels(:,1)-cut_point_base(1) acr_voxels(:,2)-cut_point_base(2) acr_voxels(:,3)-cut_point_base(3)] * SAPCA(:,3);
trim_acr_vol_id = find(sign(d) ~= sign_pos);
for id = 1:size(trim_acr_vol_id,1)
    imgd(acr_voxels(trim_acr_vol_id(id),1),acr_voxels(trim_acr_vol_id(id),2),acr_voxels(trim_acr_vol_id(id),3)) = TA_ID;
end

% trim the acromion with acr_trim_rate (to PCA)
d = [acr_voxels(:,1)-cut_point(1) acr_voxels(:,2)-cut_point(2) acr_voxels(:,3)-cut_point(3)] * SAPCA(:,3);
acr_PCA_vol_id = find(sign(d) ~= sign_pos);

% rearrange acr_PCA_vol_id so that trim_acr_vol_id takes the front
% positions
diff_id = setdiff(acr_PCA_vol_id,trim_acr_vol_id);
acr_PCA_vol_id = [trim_acr_vol_id; diff_id];

fprintf('done in %4.4f s\n',toc);

% step 2 form conoid
mean_acr = sum(acr_voxels(acr_PCA_vol_id,1:3)) / size(acr_PCA_vol_id,1);

% 2nd cutting plane
[TAPCA,sc] = princomp(acr_voxels(acr_PCA_vol_id,1:3));

cut_point_z = median(sc(1:size(trim_acr_vol_id,1),PCA_axis));
cut_point = mean_acr + TAPCA(:,PCA_axis)'*cut_point_z;

% debug
% cut_point
if (render_con_output == 1)
cut_point_cube = [cut_point; cut_point+[0 1 0]; cut_point+[1 0 0]; cut_point+[0 -1 0]; cut_point+[-1 0 0]];
figure(10)
hold on
[t]=delaunay(cut_point_cube(:,1),cut_point_cube(:,2));
trisurf(t,cut_point_cube(:,1),cut_point_cube(:,2),cut_point_cube(:,3),'facecolor',[0 0 1],'edgecolor',[0.8 0.8 0.8]);
end

minx = max(min([apex(1); acr_voxels(trim_acr_vol_id,1)]) - 3,1);
miny = max(min([apex(2); acr_voxels(trim_acr_vol_id,2)]) - 3,1);
minz = max(min([apex(3); acr_voxels(trim_acr_vol_id,3)]) - 3,1);
maxx = min(max([apex(1); acr_voxels(trim_acr_vol_id,1)]) + 3,s(1));
maxy = min(max([apex(2); acr_voxels(trim_acr_vol_id,2)]) + 3,s(2));
maxz = min(max([apex(3); acr_voxels(trim_acr_vol_id,3)]) + 3,s(3));

%% Conoid analysis
if (run_conoid_analysis == 1)
con_output0 = [];

for setting = 1:2
    if (setting == 1)
        conoid_sv_precision = max_conoid_sv_precision;
    elseif (setting == 2)
        conoid_sv_precision = 2;
    end
    
    tic,
    fprintf('Segmenting conoid...');
    
    % convert to subvoxel metric space
    con_bound_sv = zeros((maxx-minx)*conoid_sv_precision,(maxy-miny)*conoid_sv_precision,(maxz-minz)*conoid_sv_precision,'uint8');
    apex_sv = (apex-[minx miny minz]) * conoid_sv_precision;
    lut1 = double(1:(maxx-minx)*conoid_sv_precision);
    lut1 = round(lut1/conoid_sv_precision + minx);
    lut2 = double(1:(maxy-miny)*conoid_sv_precision);
    lut2 = round(lut2/conoid_sv_precision + miny);
    lut3 = double(1:(maxz-minz)*conoid_sv_precision);
    lut3 = round(lut3/conoid_sv_precision + minz);

    con_bound_sv(:,:,:) = imgd(lut1,lut2,lut3);
    % convert trimmed acr volume to subvoxel metric space
    trim_acr_sv = zeros(size(trim_acr_vol_id,1) * (conoid_sv_precision^3),3);
    id = 0;
    for x=1:size(con_bound_sv,1)
        for y=1:size(con_bound_sv,2)
            for z=1:size(con_bound_sv,3)
                v = con_bound_sv(x,y,z);
                if (v == TA_ID)
                    id = id + 1;
                    trim_acr_sv(id,:) = [x y z];
                end
            end
        end
    end
    trim_acr_sv = trim_acr_sv(1:id,:);
    clear lut1 lut2 lut3;

    cut_point_sv = (cut_point - [minx miny minz]) * conoid_sv_precision;
    sc_sv = [trim_acr_sv(:,1)-cut_point_sv(1) trim_acr_sv(:,2)-cut_point_sv(2) trim_acr_sv(:,3)-cut_point_sv(3)]*TAPCA;

    trim_acr_med_plane_sv_id = find(PRECISION >= sc_sv(:,PCA_axis) & sc_sv(:,PCA_axis) >= -PRECISION);
    %ddd = sc_sv(trim_acr_med_plane_sv_id,setdiff(1:3,PCA_axis));
    zprj = round(sc_sv(trim_acr_med_plane_sv_id,setdiff(1:3,PCA_axis)));
    clear sc_sv;
    zprj_u_loc = unique(zprj, 'rows');
    conoid_base_area = size(zprj_u_loc,1);
    %[K,conoid_base_area] = convhull(zprj_u_loc(:,1),zprj_u_loc(:,2));
    conoid_height = abs([apex_sv(1)-cut_point_sv(1) apex_sv(2)-cut_point_sv(2) apex_sv(3)-cut_point_sv(3)] * TAPCA(:,PCA_axis));
    con_voxels_sv = zeros(conoid_base_area * conoid_height, 3,'int16');
    size(con_voxels_sv)
    conoid_base_area = conoid_base_area / ((cvt_precision * conoid_sv_precision)^2);
    conoid_height = conoid_height / (cvt_precision * conoid_sv_precision);

    % conoid
    id_c = 0;
    id_c2 = 0;
    cache_size = uint16(size(trim_acr_med_plane_sv_id,1) / 16 + 1);
    cache = zeros(cache_size * conoid_height * 4, 3,'int16');
    line_precision = 0;

    for id = 1:size(trim_acr_med_plane_sv_id,1)
        x = trim_acr_sv(trim_acr_med_plane_sv_id(id),1);
        y = trim_acr_sv(trim_acr_med_plane_sv_id(id),2);
        z = trim_acr_sv(trim_acr_med_plane_sv_id(id),3);

        [xs ys zs] = bresenham_line3d(apex_sv,[x y z],line_precision);
        add = unique(round([xs' ys' zs']), 'rows');

        if (size(add,1) > 1)
            cache((id_c+1):(id_c+size(add,1)),:) = add;
            id_c = id_c + size(add,1);
        end

        % to avoid "out of mem" issue
        if (mod(id,cache_size)==0)
            fprintf('.');
            cache2 = unique(cache(1:id_c,:),'rows');
            con_voxels_sv((id_c2+1):(id_c2+size(cache2,1)),:) = cache2;
            id_c2 = id_c2 + size(cache2,1);
            id_c = 0;
        end
    end
        if (mod(id,cache_size)~=0)
            fprintf('.');
            cache2 = unique(cache(1:id_c,:),'rows');
            con_voxels_sv((id_c2+1):(id_c2+size(cache2,1)),:) = cache2;
            id_c2 = id_c2 + size(cache2,1);
            id_c = 0;
        end
    clear cache cache2 trim_acr_sv;
    con_voxels_sv = con_voxels_sv(1:id_c2,:);
    size(con_voxels_sv)
    con_voxels_sv = unique(con_voxels_sv,'rows');

    % segment conoid
    trim_acr_vol = size(trim_acr_vol_id,1) / (cvt_precision^3);
    sub_acr_vol = 0;
    conoid_vol = 0;
    conoid_acr_vol = 0;
    conoid_hum_vol = 0;

    for id = 1:size(con_voxels_sv,1)
        x_sv = con_voxels_sv(id,1);
        y_sv = con_voxels_sv(id,2);
        z_sv = con_voxels_sv(id,3);

        con_bound_sv(x_sv,y_sv,z_sv) = con_bound_sv(x_sv,y_sv,z_sv) + SA_ID;

        if (con_bound_sv(x_sv,y_sv,z_sv) == SA_ID)
            sub_acr_vol = sub_acr_vol+1;
        end

        if (con_bound_sv(x_sv,y_sv,z_sv) == CA_ID)
            conoid_acr_vol = conoid_acr_vol+1;
        end

        if (con_bound_sv(x_sv,y_sv,z_sv) == CH_ID)
            conoid_hum_vol = conoid_hum_vol+1;
        end
        conoid_vol = conoid_vol+1;
    end

    % correct boundary error
    con_surf_voxels_sv = zeros(size(con_voxels_sv,1),3);
    id_cs = 0;
    id_sa = 0;
    v_cs_surf_sv = 0;
    v_sa_surf_sv = 0;
    v_ca_surf_sv = 0;
    v_ch_surf_sv = 0;
    for id = 1:size(con_voxels_sv,1);
        x = con_voxels_sv(id,1);
        y = con_voxels_sv(id,2);
        z = con_voxels_sv(id,3);
        if (((con_bound_sv(x,y,z-1) < SA_ID) || (con_bound_sv(x,y,z+1) < SA_ID) || (con_bound_sv(x,y-1,z) < SA_ID) || (con_bound_sv(x,y+1,z) < SA_ID) || (con_bound_sv(x-1,y,z) < SA_ID) || (con_bound_sv(x+1,y,z) < SA_ID)))
            id_cs = id_cs + 1;
            con_surf_voxels_sv(id_cs,:) = [x y z];
            v = (sum(sum(sum(con_bound_sv(x-1:x+1,y-1:y+1,z-1:z+1) < SA_ID))))/27;
            v_cs_surf_sv = v_cs_surf_sv + max([min([v 1]) 0]);
            if (con_bound_sv(x,y,z) == SA_ID)
                v_sa_surf_sv = v_sa_surf_sv + v;
            end

            if (con_bound_sv(x,y,z) == CA_ID)
                id_sa = id_sa + 1;
                v_ca_surf_sv = v_ca_surf_sv + v;
            end

            if (con_bound_sv(x,y,z) == CH_ID)
                id_sa = id_sa + 1;
                v_ch_surf_sv = v_ch_surf_sv + v;
            end
        end
    end
    con_surf_voxels_sv = con_surf_voxels_sv(1:id_cs,:);
    conoid_vol = (conoid_vol - v_cs_surf_sv - size(trim_acr_med_plane_sv_id,1)*0.5) / ((cvt_precision * conoid_sv_precision)^3);
    sub_acr_vol = (sub_acr_vol - v_sa_surf_sv) / ((cvt_precision * conoid_sv_precision)^3);
    conoid_acr_vol = (conoid_acr_vol - size(trim_acr_med_plane_sv_id,1)*0.5) / ((cvt_precision * conoid_sv_precision)^3);
    conoid_hum_vol = (conoid_hum_vol - v_ch_surf_sv) / ((cvt_precision * conoid_sv_precision)^3);
    clear con_bound_sv hum_voxels acr_voxels con_voxels_sv;

    elapsed_time = toc;
    fprintf('done in %4.4f s\n',elapsed_time);
    acc = ((conoid_vol*3/conoid_base_area - conoid_height)/conoid_height*100);
    con_output0 = [con_output0; conoid_base_area conoid_height trim_acr_vol conoid_vol sub_acr_vol conoid_acr_vol conoid_hum_vol acc elapsed_time];
    
end
end

%% MAHD analysis
if (run_MAHD_analysis == 1)
    tic,
    fprintf('Finding MAHD ...');
    MAHD_output0 = [];

    MAHD_sv_precision = max_MAHD_sv_precision;
    
    % convert to subvoxel metric space
    con_bound_sv = zeros((maxx-minx)*MAHD_sv_precision,(maxy-miny)*MAHD_sv_precision,(maxz-minz)*MAHD_sv_precision,'uint8');
    apex_sv = (apex-[minx miny minz]) * MAHD_sv_precision;
    lut1 = double(1:(maxx-minx)*MAHD_sv_precision);
    lut1 = round((lut1-1)/MAHD_sv_precision + minx + 0.5);
    lut2 = double(1:(maxy-miny)*MAHD_sv_precision);
    lut2 = round((lut2-1)/MAHD_sv_precision + miny + 0.5);
    lut3 = double(1:(maxz-minz)*MAHD_sv_precision);
    lut3 = round((lut3-1)/MAHD_sv_precision + minz + 0.5);

    con_bound_sv(:,:,:) = imgd(lut1,lut2,lut3);
    clear lut1 lut2 lut3;
    
    n_voxel_sv = size(con_bound_sv,1) * size(con_bound_sv,2) * size(con_bound_sv,3);
    hum_voxels_sv = zeros(n_voxel_sv/27,3);
    acr_voxels_sv = zeros(n_voxel_sv/27,3);

    id_h = 0;
    id_a = 0;
    for x=2:step:size(con_bound_sv,1)-1;
        for y=2:step:size(con_bound_sv,2)-1;
            for z=2:step:size(con_bound_sv,3)-1;
                v = con_bound_sv(x,y,z);
                if (v ~= 0)
                    if (((v ~= con_bound_sv(x,y,z-1)) || (v ~= con_bound_sv(x,y,z+1)) || (v ~= con_bound_sv(x,y-1,z)) || (v ~= con_bound_sv(x,y+1,z)) || (v ~= con_bound_sv(x-1,y,z)) || (v ~= con_bound_sv(x+1,y,z))))
                        if (v == 2)
                            id_h = id_h + 1;
                            hum_voxels_sv(id_h,:) = [x y z];
                        else if (v == TA_ID)
                            id_a = id_a + 1;
                            acr_voxels_sv(id_a,:) = [x y z];
                            end
                        end
                    end
                end
            end
        end
    end
    hum_voxels_sv = hum_voxels_sv(1:id_h,:);
    acr_voxels_sv = acr_voxels_sv(1:id_a,:);
    
    min_d = size(con_bound_sv,1)^2 + size(con_bound_sv,2)^2 + size(con_bound_sv,3)^2;
    idA = -1;
    for i1 = 1:size(hum_voxels_sv,1)
        for i2 = 1:size(acr_voxels_sv,1)
            d = (hum_voxels_sv(i1,1) - acr_voxels_sv(i2,1))^2 + (hum_voxels_sv(i1,2) - acr_voxels_sv(i2,2))^2 + (hum_voxels_sv(i1,3) - acr_voxels_sv(i2,3))^2;
            if (d < min_d)
                min_d = d;
                idA = i2;
                idH = i1;
            end
        end
    end
    min_d = min_d^0.5 / (cvt_precision * MAHD_sv_precision);

    if (idA ~= -1)
        A = ((acr_voxels_sv(idA, :)-1) / MAHD_sv_precision + [minx+1 miny+1 minz+1]);
        H = ((hum_voxels_sv(idH, :)-1) / MAHD_sv_precision + [minx+1 miny+1 minz+1]);
    else
        error('\n MAHD analysis failed!\n')
    end
    
    if (run_MAHD_analysis == 1)
        [x y z] = bresenham_line3d(A,H,2);
        figure(10)
        hold on
        plot3(x,y,z);
    end
    
    elapsed_time = toc;
    fprintf('done in %4.4f s\n',elapsed_time);
    clear con_bound_sv hum_voxels_sv acr_voxels_sv;
    A = round((A-1) ./ (scale_factor * cvt_precision) + 0.99);
    H = round((H-1) ./ (scale_factor * cvt_precision) + 0.99);
    
    % convert to coordinate system of drawing tool
    A = [os(1)-A(1) os(2)-A(2) A(3)-1];
    H = [os(1)-H(1) os(2)-H(2) H(3)-1];
        
    MAHD_output0 = [MAHD_output0; A(1) A(2) A(3) H(1) H(2) H(3) min_d elapsed_time];
end

%% Output
    i = 1;
    fprintf('\nOutputs\n');
    if (run_conoid_analysis == 1)
        con_output = [con_output; con_output0(i,:)];
        
        fprintf('\tconoid_base_area\t%.2f\tmm2\n',con_output0(i,1));
        fprintf('\tconoid_height\t%.2f\tmm\n',con_output0(i,2));
        fprintf('\ttrim_acr_vol\t%.2f\tmm3\n',con_output0(i,3));
        fprintf('\tconoid_vol\t%.2f\tmm3\n',con_output0(i,4));
        fprintf('\tsub_acr_vol\t%.2f\tmm3\n',con_output0(i,5));
        fprintf('\tconoid_acr_vol\t%.2f\tmm3\n',con_output0(i,6));
        fprintf('\tconoid_hum_vol\t%.2f\tmm3\n',con_output0(i,7));
        fprintf('\tsegmentation_error\t%.2f%%\t\n',con_output0(i,8));
        fprintf('\tsegmentation_time\t%.2fs\t\n',con_output0(i,9));
    end
    if (run_MAHD_analysis == 1)
        MAHD_output = [MAHD_output; MAHD_output0(i,:)];
        
        fprintf('\tacr_x\t%.2f\t\n',MAHD_output0(i,1));
        fprintf('\tacr_y\t%.2f\t\n',MAHD_output0(i,2));
        fprintf('\tacr_z\t%.2f\t\n',MAHD_output0(i,3));
        fprintf('\thum_x\t%.2f\t\n',MAHD_output0(i,4));
        fprintf('\thum_y\t%.2f\t\n',MAHD_output0(i,5));
        fprintf('\thum_z\t%.2f\t\n',MAHD_output0(i,6));
        fprintf('\tMAHD\t%.2f\tmm\n',MAHD_output0(i,7));
        fprintf('\tsegmentation_time\t%.2fs\t\n',MAHD_output0(i,8));
    end
    if (run_sca_ver_analysis == 1)
        sca_ver_output = [sca_ver_output; [S2V_angles_deg R(1,:) R(2,:) R(3,:)]];
        fprintf('\ty rotation amount: \t%.2fdegrees\t\n',S2V_angles_deg(1));
        fprintf('\tx rotation amount: \t%.2fdegrees\t\n',S2V_angles_deg(2));
        fprintf('\tz rotation amount: \t%.2fdegrees\t\n',S2V_angles_deg(3));
        fprintf('\trotation_matrix\n');
        R,
    end
    
    % render outputs
    tic
    fprintf('\nRendering surfaces...\n');
    
    figure(8)
    axis equal;
    view(3);
    
    figure(10);
    axis equal;
    view(3);
    
    if (render_con_output == 1 && run_conoid_analysis == 1)
        % conoid
        try
            p = unique([con_surf_voxels_sv(:,1) / conoid_sv_precision + minx con_surf_voxels_sv(:,2) / conoid_sv_precision + miny con_surf_voxels_sv(:,3) / conoid_sv_precision + minz],'rows');
            if (size(p) < 32)
                throw;
            else
                [t]=MyCrustOpen(p);
            end

            figure(8)
            hold on
            trisurf(t,p(:,1),p(:,2),p(:,3),'facealpha',0.25,'facecolor',yellow,'edgecolor',[0.8 0.8 0.8]);
        catch exception
            fprintf('\nRendering conoid surface failed!\n');
        end
    end
    
    if ((render_con_output == 1 && run_conoid_analysis == 1) || run_MAHD_analysis == 1)
        max_complexity = 1024;
        for id = [AC_ID HU_ID SC_ID VE_ID TA_ID]
            try
                ids = find(surf_voxels(:,4) == id);
                step = max(1,round(size(ids) / max_complexity));
                p = surf_voxels(ids(1:step:size(ids)),1:3);
                if (size(p) < 32)
                    throw;
                else
                    [t]=MyCrustOpen(p);
                end

                figure(10)
                hold on
                fprintf('\tGenerating surface of label id %d...\n',id);
                trisurf(t,p(:,1),p(:,2),p(:,3),'facealpha',0.25,'facecolor',facecolor(id,:),'edgecolor',[0.8 0.8 0.8]);
            catch exception
                fprintf('\nRendering surface of label id %d failed!\n',id);
            end
        end
        for id = [TA_ID] % figure 8 is dedicated for showing a particular bone
            try
                ids = find(surf_voxels(:,4) == id);
                step = max(1,round(size(ids) / max_complexity));
                p = surf_voxels(ids(1:step:size(ids)),1:3);
                if (size(p) < 32)
                    throw;
                else
                    [t]=MyCrustOpen(p);
                end

                figure(8)
                hold on
                fprintf('\tGenerating surface of label id %d...\n',id);
                trisurf(t,p(:,1),p(:,2),p(:,3),'facealpha',0.25,'facecolor',facecolor(id,:),'edgecolor',[0.8 0.8 0.8]);
            catch exception
                fprintf('\nRendering surface of label id %d failed!\n',id);
            end
        end
        h = figure(10);
        saveas(h,[filename '_10.fig']);
        h = figure(8);
        saveas(h,[filename '_8.fig']);
    end
    if (run_sca_ver_analysis == 1)
        h = figure(9);
        saveas(h,[filename '_9.fig']);
    end
    fprintf('done in %.2f s\n',toc);
end

% print batch conoid analysis outputs
if (run_conoid_analysis == 1)
    fprintf('\nConoid Analysis Outputs\n');
    for i=1:size(con_output,1)
       fprintf('%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\n',allfiles(i).name,con_output(i,:)); 
    end
end

% print batch MAHD analysis outputs
if (run_MAHD_analysis == 1)
    fprintf('\nMAHD Analysis Outputs\n');
    for i=1:size(MAHD_output,1)
       fprintf('%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\n',allfiles(i).name,MAHD_output(i,:)); 
    end
end

% print batch Scapula - Vertebra analysis outputs
if (run_sca_ver_analysis == 1)
    fprintf('\nScapula Vertebra Analysis Outputs\n');
    for i=1:size(sca_ver_output,1)
       fprintf('%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\n',allfiles(i).name,sca_ver_output(i,:)); 
    end
end