% Plot PCA axis

% Origin of axis is at mean surface crd
% Axis are scaled suitably using min/max crd

% -------------------------------------------------------------------------

function [fighndl]=plot3axes(x,evector,name,fighndl)

% Input
% x(nv,3)               Vertex coordinates of surface
% evector(3,3)          Eigenvector from PCA : col are axis
% name                  Name of image
% fighndl               Figure handle of existing figure

% Output
% fighndl               Figure handle 

% Local

% -------------------------------------------------------------------------

if (~isempty(evector))
    % Display axis
    value=exist('fighndl','var');
    if ( value==0 | isempty(fighndl))
        fighndl=figure('Name',name);  % new figure
        title(name);
        legend('Primary','Secondary','Tertiary',4);   % 4 = bottom right  
    else
        figure(fighndl);  % existing figure
%         legend('1','2','3',4);   % 4 = bottom right  % commented as superimposing plots causes wrong legend
        hold on
    end
    
    xmean=mean(x);
    xmin=min(x);
    xmax=max(x);
    sf=norm(xmax-xmin);  % diagonal of bounding cube
    sf=0.5*sf;
    
    x1=[xmean(1);xmean(1)+sf*evector(1,1)];
    y1=[xmean(2);xmean(2)+sf*evector(2,1)];
    z1=[xmean(3);xmean(3)+sf*evector(3,1)];
    x2=[xmean(1);xmean(1)+sf*evector(1,2)];
    y2=[xmean(2);xmean(2)+sf*evector(2,2)];
    z2=[xmean(3);xmean(3)+sf*evector(3,2)];
    x3=[xmean(1);xmean(1)+sf*evector(1,3)];
    y3=[xmean(2);xmean(2)+sf*evector(2,3)];
    z3=[xmean(3);xmean(3)+sf*evector(3,3)];

    plot3(x1,y1,z1,'r-',x2,y2,z2,'g-',x3,y3,z3,'k-')  % red, green, black
        view(45,30) 
        axis tight 
        daspect([1,1,1])

    hold off
else
    % No axis : do not plot
    fighndl=[];
end

% -------------------------------------------------------------------------
% End of Function