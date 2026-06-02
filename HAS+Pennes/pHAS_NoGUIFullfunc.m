% script used to run HAS without using the user interface.
% see HAS_NoGUI_Instructions.txt for instructions
% Sara Johnson, February 2020
% with modifications by Michelle Kline, April 2020
%   - call CalcHAS_ERFA8e instead of duplicating calculations here
% 
% Michelle Kline, July 2020
% added required variables for phase correction, if usePhaseCorrection == true or
% usePhaseTimeRev == true

%%% VARIABLES THAT MUST BE SET BEFORE CALLING THIS SCRIPT %%%
% ----- pathname and filename of ERFA and model files -----
% erfaFullFile   = full pathname/filename of ERFA .mat file
% modelFullFile  = full pathname/filename of model .mat file
% ----- variables used for the HAS calculations -----
% c0        = Speed of sound in water (m/s)
% rho0      = Density of water 
% Dx        = Modl resolution in 2nd dimension
% Dy        = Modl resolution in 1st dimension
% Dz        = Modl resolution in 3rd dimension
% c         = vector of Speed of sound values in each Modl domain (m/s)
% a         = vector of total Attenuation in each Modl domain 
% rho       = vector of Density in each Modl domain 
% Pr        = Acoustic power output of transducer
% Modl      = 3D segmented model, corresponds to values in a, rho, c
% perfa     = Input the perfa separately because it takes a long time to load
% ----- for the following inputs, a positive (+) value progresses to a higher slice number -----
% offsetxmm  = mechanical offset from center of Modl (along 2nd dimension) (mm)
% offsetymm  = mechanical offset from center of Modl (along 1st dimension) (mm)
% dmm        = Distance from Xducer base to Modl base (mm)
% hmm        = electronic steering in x-direction (along 2nd dimension) (mm)
% vmm        = electronic steering in y-direction (alond 1st dimension) (mm)
% zmm        = electronic steering in z-direction (along 3rd dimension) (mm)
% ----- options that would be set by user in the GUI -----
% usePhaseCorrection = use phase correction? (true or false)
% usePhaseTimeRev    = use phases from time reversal? (true or false)
% useRayTracing      = use ray tracing? (true or false)
% numrefl    = number of reflections desired
% ----- pathname and filename of output .mat files -----
% poutFullFile  
% QFullFile  
% ppFullFile  

%%% OUTPUTS %%%
% Q         = SAR (power deposition) in 3D Modl
% maxQ      = maximum Q-value in 3D Modl
% pout      = pressure pattern in 3D Modl
% pp        = pressure on input plane

%%% LOAD ERFA, MODEL, and (OPTIONAL) PHASE CORRECTION FILES %%%
% if ERFA and model files are already loaded, don't attempt to re-load
function [pout,Z,angpgvect,phasematsv,pR,pHAS] = pHAS_NoGUIFullfunc(Modl,pHAS,pR,pERFA)

if ~isfield(pR,'ERFA_load')
    pR.ERFA_load = 1;
end
if ~isfield(pR,'modl_load')
    pR.modl_load = 1;
end
if pR.ERFA_load == 1
    disp('Loading ERFA and model files...');
    if ~isfield(pR,'erfaFullFile')
        disp('Please set erfaFullFile (ERFA file location)');
        return;
    else
        load(pR.erfaFullFile);
    end
    disp('Finished loading ERFA.');
end
if pR.modl_load == 1
    if ~isfield(pR,'modelFullFile')
        disp('Please set modelFullFile (model file location)');
        return;
    else
        load(pR.modelFullFile);
    end
    disp('Finished loading model.');
end
if pR.usePhaseCorrection == true || pR.usePhaseTimeRev == true
    if ~exist('phase_load','var')
        phase_load = 1;
    end
    if phase_load == 1
        if ~isfield(pR,'phaseFullFile')
            disp('Please set phaseFullFile (phase correction file location)');
            return;
        else
            %load(pR.phaseFullFile);
        end
        %disp('Finished loading phase corrections.');
    end
end
pR.AmpPh1 = pR.Phases(pR.jXducerSection).AmpPh1;
if pR.wbr; disp('Finished loading phase corrections.'); end
%%% INPUT VARIABLE UNIT CONVERSION %%%
pR.h=pHAS.hmm/1000; % convert to m units.
pR.v=pHAS.vmm/1000; 
pR.z=pHAS.zmm/1000; 
pR.d=pHAS.dmm/1000;
pR.f=pERFA.fMHz*1e6;	% convert to Hz.

%%% ERFA PARAMETERS %%%
% perfa = pERFA.perfa;
% if strcmp(class(perfa),'double')  %#ok<STISA>
%     perfa=single(perfa); 
% end % always use single precision perfa to save memory.
% 
% pR.Rmm=pERFA.R*1000; % R in m for most calcs, Rmm in mm for gui display.
% %pR.Rmm = Rmm;
% [lmaxerfa,mmaxerfa,~]=size(perfa);  % size of ERFA plane.
% dyerfa=pERFA.Len(1)/(lmaxerfa-1); 
% dxerfa=pERFA.Len(2)/(mmaxerfa-1); % sample spacing in ERFA plane, in m.
% Dyerfa=dyerfa*1000; 
% Dxerfa=dxerfa*1000;    
% 
% yaxiserfa=Dyerfa*(-(lmaxerfa-1)/2:(lmaxerfa-1)/2);    % setting up axes in mm units for interpolation.
% xaxiserfa=Dxerfa*(-(mmaxerfa-1)/2:(mmaxerfa-1)/2); 
% xaxiserfaoffs=xaxiserfa+pHAS.offsetxmm; 
% yaxiserfaoffs=yaxiserfa+pHAS.offsetymm;   % adjust axes for offsets, all in mm.
% 
%%% MODEL PARAMETERS %%%%%%%%%%
% Copied from loadmodel6.m
% sm=size(Modl);
% % Copied from drawmodel6.m
% xaxisinterp=Dx*(-(sm(2)-1)/2:(sm(2)-1)/2); % axes (between centers of end points) for imagesc and interp.
% yaxisinterp=Dy*(-(sm(1)-1)/2:(sm(1)-1)/2);     % Dx, Dy and Dz all in mm.
% %zaxis=Dz*(1:sm(3)); % longitudinal axis has full Dz at the center of the first voxel, since HAS calculates
%    % travel through a full distance Dz for each voxel and attributes the resulting pressure to that voxel.
% %xaxis=xaxisinterp; yaxis=yaxisinterp;   % to allow use of legacy xaxis and y axis labels.
% Lx=(sm(2)-1)*Dx; 
% Ly=(sm(1)-1)*Dy; 
% Lz=sm(3)*Dz; 
% LenModl=[Lx,Ly,Lz]; % in mm. Note: Lz is overall length.
% lx=Lx/1000; 
% ly=Ly/1000; % convert to m units.
% Dv=[Dx Dy Dz];  %  Model resolution (dim 1, 2, 3)) IN PARAMETER FILE
% 
% [lmax,mmax,nmax] = size(Modl);  % The size of the model (lmax,mmax,nmax) sets the size of the simulation
%     % space.  lmax and y are vertical; mmax and x are horizontal.  lmax and mmax are therefore also the size of
	% the pressure pattern pp on the front plane of the Modl, interpolated from perfa. Note: lmax and mmax 
	% should be ODD numbers to keep fftshift symmetrical with the dc term at exact center of spectrum.
    % Note that lmax,mmax,nmax are duplicates of sm(1),sm(2),sm(3) because of legacy. 

% MMK FIXME - ask Sara for comment on aabs
% if ~exist('aabs','var')
%     aabs=   pHAS.a_abs;     %a; 
% end
   
% run the HAS calculations
if pR.wbr; disp('Running simulation...'); end
%usingGUI = false;
%***********************************************************************
%CalcHAS_ERFA8eFullnoGUI;         %THIS IS THE HAS KERNEL !!! 
%DLP 1/31/22  set this 
%back to use the correct radial distance "rave"
[pout,Z,angpgvect,phasematsv,pR,pHAS] = CalcHAS_ERFA8eFullnoGUIfunc(Modl,pHAS,pR,pERFA);
%***********************************************************************

if pR.wbr; disp('Simulation complete'); end

HAS_savefile = pR.HAS_savefile; % zach j addded
if ~exist('HAS_savefile','var')
    HAS_savefile = 1;
end
if HAS_savefile == 1
    % save outputs to .mat files
    if pR.wbr; disp('Saving output files...'); end
    %Savefiles7
    if pR.wbr; disp('Saving output complete'); end
end