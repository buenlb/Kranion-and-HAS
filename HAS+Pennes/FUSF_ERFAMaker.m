%  This is a modified version of D.A. Christiensen's ERFA maker. 
%  Creates 7 efra planes, one for each section of
%  a hemispherical trasducer. Geometry loaded from Kranion. 
%
% @INPUTS
%   Seg: A stucture containing transducer element locations divied into 7
%   sections. 
%   TH_rotation: 7 rotation angles about the X axis to bring each of the 7
%   transducer segments to align with the axis of propagation. 
%   PHI_rotation: 7 rotation angles about the Y axis to bring each of the
%   7 transducer segments to align with the axis of propagation. 
%   fullarray: A maxrix of all trasducer element locations.
%   fullfilepath: Path to the Kranion export (.krx) file in use. 
%
%   Zach Johnson and Taylor Webb 
%   June 5, 2024
    
function FUSF_ERFAMaker(Params,Seg,TH_rotation,PHI_rotation,fullarray,fullfilepath)
disp('No precomputed ERFA file found...computing ERFA now.')

% Set ERFA Parameters from Params struct in the variable space
fields = fieldnames(Params);
for i = 1:length(fields)
    eval([fields{i} ' = Params.' fields{i} ';']);
end

Dim=[Dv,Dh]; anginc=[imax,kmax]; Len=[Lv,Lh]; planinc=[lmax,mmax]; % combine some params.

%% Launch ERFA Maker for each array section
%Run ERFA for each segment, rotating each by:
% [TH_centers_rot, PHI_centers_rot]

if 2*round(lmax/2)==lmax || 2*round(mmax/2)==mmax
  hdb=warndlg('The number of ERFA plane increments should be ODD integers. Continue only if even is okay.',...
     'Not Odd', 'modal'); uiwait(hdb); 
end

f=fMHz*1e6;	% convert to Hz.

i=1:imax; k=1:kmax; l=1:lmax; m=1:mmax;     % set up indices.  Careful: i used as index here, j is imag number.
dth=single(Dv/(R*imax)); dphi=single(Dh/(R*kmax));		% incremental size of source angle (in radians).
dyp=single(Lv/(lmax-1)); dxp=single(Lh/(mmax-1));	% incremental size of steps in ERFA plane (in m).
th1=dth*(i-round(imax/2));      % angle row vector, centered; imax and kmax should be odd for symmetry.
phi1=dphi*(k-round(kmax/2));	% angle row vector, centered.
thmesh1=repmat(th1',1,kmax); % imax x kmax matrix of theta values, 'meshgrid' style, over entire rectangular
                                % angular area that encompasses xducer.
phimesh1=repmat(phi1,imax,1);% imax x kmax matrix of phi values, 'meshgrid' style, over entire xducer area.

Zm=rho0*c0;     % impedance of medium that waves radiate into; rho0 and c0 read in the parameter file.


for ns = 1:length(Seg)
    % Rotate array segment into correct orientation, then save spherical coords
    % into thvect and phivect.

    Ry =        [cos(PHI_rotation(ns)) 0 sin(PHI_rotation(ns)); 
                0 1 0; 
                -sin(PHI_rotation(ns)) 0 cos(PHI_rotation(ns));];

% Create the rotation matrix about X' axis by psi
     Rx =        [1 0 0;
                 0 cos(TH_rotation(ns)) sin(TH_rotation(ns));
                  0 -sin(TH_rotation(ns)) cos(TH_rotation(ns));];
   
    a = Seg(ns).ElemLoc_Cart(:,1);
    b = Seg(ns).ElemLoc_Cart(:,2);
    c = Seg(ns).ElemLoc_Cart(:,3);
    rotated = [Rx*Ry*[a b c]']';
    [thvect, phivect] = cart2sphd(rotated(:,1),rotated(:,2),rotated(:,3), R);
    h = figure; h.Position = [1252          87         667         890];
    %plotTx(a,b,c,h)
    %hold on
    %plotTx(rotated(:,1),rotated(:,2),rotated(:,3),h);
    %view([1,0,0]);
    
    close(h);
    % a = Seg(ns).ElemLoc_Cart(:,1);
    % b = Seg(ns).ElemLoc_Cart(:,2);
    % c = Seg(ns).ElemLoc_Cart(:,3);
    % rotated = [Ry*Rz*[a b c]']';
    % [phivect, thvect,~] = cart2sph(rotated(:,1),rotated(:,2),rotated(:,3));
    
        close all
        plot3(fullarray(:,1),fullarray(:,2),fullarray(:,3),'g*');
        axis image, grid on, hold on, xlabel('x'),ylabel('y'),zlabel('z')
        plot3(a,b,c,'ko');
        plot3(rotated(:,1),rotated(:,2),rotated(:,3),'r*');
    h1=waitbar(0,'Evaluating segment ERFA pressure pattern...');
    numelem=size(thvect,1);
    perfa=zeros(lmax,mmax,numelem,'single');     % preallocation for speed.         
    pointmap=zeros(imax,kmax,'single'); % initialize map of valid points inside xducer elements.


    % --- Loop through the elements (of the current segment) starting here ---
    for g=1:numelem       % cycle through elements, finding a page of ERFA for each element.
        
        distfromelem=R*acos((cos(thvect(g))*cos(thmesh1)).*cos(phivect(g)-phimesh1)...
            +sin(thvect(g))*sin(thmesh1));   % great circle distance from center of element g.

        indelemlin= find(distfromelem<=relem);  % find linear indices inside round element g.
        [indelemi, indelemk]=ind2sub([imax,kmax],indelemlin);   % convert to subscripts.
        is=min(indelemi);  ie=max(indelemi);    % find start and end indices in area around element g.
        ks=min(indelemk);  ke=max(indelemk);
        th2=th1(is:ie);     % angle row vector only around element.
        phi2=phi1(ks:ke);     % angle row vector only around element.
        isize=size(th2,2);      % size of rectangle enclosing element g only, so only integrate over element.
        ksize=size(phi2,2);
%                 if isize==1 || isize==0 || ksize==1 || ksize==0; 
%                     errordlg(['Only 0 or 1 sample points inside at least one element. Increase number'...
%                         ' of angle increments to transducer.']); return; end
        thmesh2=repmat(th2',1,ksize); % isize x ksize meshgrid of theta values, only around element.
        phimesh2=repmat(phi2,isize,1);  % isize x ksize meshgrid of phi values, only around element.
        pt=zeros(imax,kmax,'single');        % start with blank pressure over entire field.  
        %TotElemArea=(4*relem^2)*numelem;   % total area of all elements combined (legacy).
        pt(is:ie,ks:ke)=1;      % pt = 1 only where there are valid points inside element.                      

        pt2=pt(is:ie,ks:ke);    % rectangular pressure matrix just around the element.
        cth=cos(thmesh2);
        ss5=R*(cth.*sin(phimesh2));
        s=repmat(ss5,[1,1,lmax]);		% 3D array: isize x ksize x lmax.
        xp=dxp*(m-round(mmax/2));		% vector of x points along (with vers. 8) the horizontal axis.

        aa5=R*(cth.*cos(phimesh2));
        a=repmat(aa5,[1,1,lmax]);		% 3D array: isize x ksize x lmax.
        b=R-a;  clear a     % to save memory

        tt= R*sin(thmesh2);
        t=repmat(tt,[1,1,lmax]);		% 3D array of t.
        yyp(1,1,:)=dyp*(l-round(lmax/2));	% turn y points into a 'page' vector, lmax pages long.
        yp=repmat(yyp,[isize,ksize,1]);		% make 3D array.
        term1=(t-yp).^2; clear t yp      % to save memory
        term3=(d-b).^2;  clear b         % to save memory

        ppc=pt2.*cth;					% ppc is isize x ksize matrix; cth is for spherical integration.
        pc=repmat(ppc,[1,1,lmax]);		% product of pressure and cos theta now 3D array.
        kk=2*pi*f/c0;
        ppi=zeros(ksize,lmax,mmax,'single');     % prealocation for speed.

        for mi=1:mmax
             %--- Rayleigh-Summerfeld integral done next with three implicit for loops ---
            r=sqrt(term1 + term3 + (s-xp(mi)).^2);	% r is 3D array for each value of xp.
            ppi(:,:,mi)=(f*R*R*dth*dphi/c0)*sum(pc.*exp(1j*(-kk*r + (pi/2)))./r);	
        end
        perfa(:,:,g)=shiftdim(sum(ppi));

        waitbar(g/numelem)
%                 if g==1
%                     hh1=imagesc(phi1,th1,pt); axis image; axis xy;
%                     title('Angle plot of element 1 - face view'); xlabel('phi (radians)'); ylabel('theta (radians)');
%                     figure; hh2=imagesc(xp,squeeze(yyp),squeeze(abs(perfa(:,:,g)))); axis image; axis xy;
%                     title('Space plot of ERFA pressure from element 1 - face view'); 
%                     xlabel('horizontal (m)'); ylabel('vertical (m)');
%                     button6=questdlg('Do you want to continue with Rayleigh-Sommerfeld calculations?',...
%                         'CONTINUE','Yes','Quit','Yes');
%                     if strcmp(button6,'Quit'); close all; close(h1); return; end
%                 end
        pointmap=pointmap+pt;

    end         % end of g loop.

    TotalElemArea=R*R*dth*dphi*sum(sum(pointmap.*cos(thmesh1)));     % Add up areas of all points.
    ptelem=sqrt(2*Zm*1/TotalElemArea);      % pressure corresponding to 1 W total.
    perfa=perfa*ptelem;     % multiply normalized perfa by the pressure corresponding to 1 W total.

    close all; close(h1)
    hh3=imagesc(phi1,th1,pointmap); axis image; axis xy
    title('Angle plot of all elements - face view'); xlabel('phi (radians)'); ylabel('theta (radians)');

    
    phiInElem=phimesh1(pointmap>0);    % use logical indexing to find angles to points inside element.
    thetaInElem=thmesh1(pointmap>0);
    Xplot=R*cos(thetaInElem).*sin(phiInElem);
    Yplot=R*sin(thetaInElem);
    Zplot=R - R*cos(thetaInElem).*cos(phiInElem);  % put Z in ordinary direction coming out of xducer.
    hh4=plot3(Xplot,Yplot,Zplot,'.'); axis image;  axis xy % make 3D plot of points found in elements.
    title('3D space plot of all elements - face view'); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); grid on;

    figure; hh5=imagesc(xp,squeeze(yyp),abs(sum(perfa,3))); axis image; axis xy
    title('Space plot of summed ERFA pressure from transducer - face view'); 
    xlabel('horizontal (m)'); ylabel('vertical (m)');

    perfa=fliplr(perfa);    % flip perfa so the x-axis (col) matches the direction of the x-axis in HAS.
    perfa_seg{ns} = perfa;
end % end of s (segment) loop.


sx=d;
% Save key parameters in .mat files
disp('Saving ERFA8.mat file...');
save(sName,'perfa_seg','fMHz','Len','sx','R','isPA','fullarray','dxp','dyp','relem','Seg','TH_rotation','PHI_rotation','-v7.3'); 
disp('Finished saving');
end   