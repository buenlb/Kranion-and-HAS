function [Q_rel,pressar,phmat,AmpPhar,pHAS,pR] = runInsightecHAS_func2(curfilepath,pHAS,Modl_in,pR)
%[Q_rel,pressar,q_apl,pressure]
%[Q_rel,pressar,phmat,AmpPhar,pHAS,pR] = runInsightecHAS_func3_wj(Modl_in,pHAS,pR)
%runInsightecHAS_func Summary of this function goes here
%   Detailed explanation goes here
% DLP: 12/11/20 changed to output final 1024 array of amplitudes and phases

test = 0;

if ~isfield(pR,'wbr')  % the wbr field controls the verboseness of the HAS function 0 is no visible output, 1 is show waitbars and progress updates
    pR.wbr = 0;
end



load([curfilepath filesep 'ERFA8.mat']);

% added from earlier Kranion HAS version
Foc_ijk = pR.Foc_ijk;
padVal = 0;
Modl_in = padarray(Modl_in,[ones(1,3)]*padVal,1);
Foc_ijk = Foc_ijk + padVal;


test = 0;
PhCo = 2;%pR.PhCo;
pR.PhCo=PhCo;
pR.useRayTracing = false;
if pR.PhCo==1 %no correction
    usePhaseCorrection = false;  % true/false
    usePhaseTimeRev = false;  % true/false
elseif pR.PhCo==2 %IInsightec phase correction
    usePhaseCorrection = true;  % true/false
    usePhaseTimeRev = false;  % true/false
elseif pR.PhCo==3 % time reversal phase correction
    usePhaseCorrection = false;  % true/false
    usePhaseTimeRev = true;  % true/false
end
pR.usePhaseCorrection = usePhaseCorrection;
pR.usePhaseTimeRev = usePhaseTimeRev;
pHAS.padsize = padVal;

if pR.wbr; disp(['PhCo=',num2str(pR.PhCo),' usePhaseCorrection=',num2str(usePhaseCorrection),' usePhaseTimeRev=',num2str(usePhaseTimeRev)]); end

Modl = Modl_in;

Powerv = 1; % "total Power*"(W) *assuming everthing in the ERFA maker is working 

[xdnx,xdny,xdnz] = size(Modl);
xcent = ceil(xdnx/2);
ycent = ceil(xdny/2);
zcent = ceil((150-pHAS.dmm)/pHAS.Dz);


aabs = pHAS.a_abs;


pHAS.CSF_vals.a = pHAS.a;
pHAS.CSF_vals.c = pHAS.c;
pHAS.CSF_vals.a_abs = pHAS.a_abs;
pHAS.CSF_vals.rho = pHAS.rho;

phmat = 0;

%DLP 10/20/20 put put in the Modl for rotation purposes
if(pR.PhCo == 3)
    sm = size(Modl);
    focusindv=round( (pHAS.geom(2) + pHAS.vmm)/pHAS.Dy + 0.5 + (sm(1)/2) );  % y index of new focus in padded Modlpd;
    focusindh=round( (pHAS.geom(1) + pHAS.hmm)/pHAS.Dx + 0.5 + (sm(2)/2) );  % distances in mm.
    focusindz=round( (pHAS.geom(3) + pHAS.zmm)/pHAS.Dz);
    Modl_in(focusindv,focusindh,focusindz) = pHAS.indtarget; %stick with brain properties at the focus
end

pressar = single(zeros(xdnx,xdny,xdnz,7));
AmpPhFull = single(zeros(128,2,8));

% pR.XDfolloc;
% XDnum=pR.XDfolloc(end-3:end);
for jXducerSection = 1:7
    jj=jXducerSection;

    tic

    pERFA.perfa = perfa_seg{jXducerSection};
    pERFA.fMHz = fMHz;
    pERFA.Len = Len;
    pERFA.sx = sx;
    pERFA.R = R;
    pERFA.isPA = isPA;
    pERFA.ElemLoc = Seg(jXducerSection).ElemLoc_Polar;
    pERFA.dxp = dxp;
    pERFA.dyp = dyp;
    pERFA.relem = relem; 


    tmms.loefa(jj)=toc;
    % Transducer power for sector

    if jXducerSection == 1
        Pr = Powerv/4;
    else 
        Pr = Powerv/8;
    end


    % load in each phase correction for each zone
    if usePhaseCorrection == true || usePhaseTimeRev == true
        fname = [curfilepath filesep 'ERFA8.mat'];       
    end

    pR.phaseFullFile = fname;
    pR.jXducerSection = jXducerSection;
    pR.Pr = Pr;
    % ______________________________________________
    rottm=tic;
    % if jXducerSection==1
    %     Modl=Modl_in;
    % else
    %     [Modl] = rotvolpivrecenter_00(Modl_in,pivind,dx,dy,dz,thetav(jXducerSection),psiv(jXducerSection),0,1); %%%%%%%%%%%%% SLOW %%%%%%%%%%%
    % end
    TH_rot = -rad2deg(TH_rotation);
    PHI_rot = -rad2deg(PHI_rotation);
    Modl = prepHAS_Model(Modl_in,PHI_rot(jXducerSection), -TH_rot(jXducerSection),Foc_ijk,-pHAS.xtilt);
    tmms.rotop(jj)=toc(rottm);
    % ______________________________________________
    % set power = PatPower(jX)
    %run HAS-nogui
    if(pR.PhCo == 3)
        [vmax, locmax] = max(Modl(:));
        %         Modl(locmax) = Modl(locmax) - 100;
        [focusindv,focusindh,focusindz] = ind2sub(size(Modl),locmax);           %This is why we put the target in here, as the maximum value
        focusindv=focusindv + pHAS.padsize ;  % y index of new focus in padded Modlpd;
        focusindh=focusindh + pHAS.padsize ;  % y index of new focus in padded Modlpd;
        disp(['focusindv=',num2str(focusindv),' focusindh=',num2str(focusindh),' focusindz=',num2str(focusindz)]);
    end
    pR.ERFA_load = 0;
    pR.modl_load = 0;
    pR.HAS_savefile = 0;
    pR.f = pERFA.fMHz * 1e6;

    
    pHAS.vmm = pHAS.SteeredFocus(2);
    pHAS.hmm = pHAS.SteeredFocus(1);
    pHAS.zmm = pHAS.SteeredFocus(3);

    pHAS.dmm = 150-Foc_ijk(3) -1;

    %pR.AmpPh1 = AmpPh1;
    %********************************************************
    %pHAS_NoGUIFull;         %THIS IS THE MAIN HAS PROGRAM !!!
    pnogi=tic;
    [pout,Z,angpgvect,phasematsv,pR,pHAS] = pHAS_NoGUIFullfunc(Modl,pHAS,pR,pERFA);
    tmms.pnogi(jj)=toc(pnogi);
    %********************************************************
    if jXducerSection == 1
        Z_modl=Z;
        %absmodl=pHAS.a_abs(Modl_in)*1e2*fMHz;     % DLP added 10/12/20 a(i) is pressure total attenuation coefficient (assume linear freq dep).
        %absmodl_modl=a_abs(Modl)*1e2*pERFA.fMHz;    % aabs(i) is pressure absorption coefficient (no random variation in it now).
        %absmodl_modl=absmodl;
    end

    if test ==1
        testPhCo;
    end
    tic
    %pout=rotvolpivrecenterinterp(pout,pivind,Dx,Dy,Dz,thetav(jXducerSection),psiv(jXducerSection),0,1);
    pout=rotvolpivrecenterinterp(pout,[ceil(xcent),ceil(ycent),Foc_ijk(3)],1,1,1,PHI_rot(jXducerSection),TH_rot(jXducerSection),0,1);
    if test == 1
        %pout2 = rotvolpivrecenterinterp(pout2,pivind,Dx,Dy,Dz,thetav(jXducerSection),psiv(jXducerSection),0,1);
        pout2 = rotvolpivrecenterinterp(pout2,[ceil(xcent),ceil(ycent),Foc_ijk(3)],1,1,1,PHI_rot(jXducerSection),TH_rot(jXducerSection),0,1);
    end
    tmms.rtitp(jj)=toc;
    %pout=rotvolpivrecenterinterp(pout,pivs,Dx,Dy,Dz,thetas,psis,0,1);
    %Q=rotvolpivrecenterinterp(Q,pivs,Dx,Dy,Dz,thetas,psis,0,1);   % now interpolated.
    if(jXducerSection==1)
        pressure = pout;
        if test == 1
            pressure2 =    pout2;
        end
        %DLP debug
        if(pR.PhCo ==3)
            [npx,npy,nel] = size(phasematsv);
            AmpPhFull(:,:,1) = pR.AmpPh1(1:128,:);
            AmpPhFull(:,:,2) = pR.AmpPh1(129:256,:);
            phmat = single(zeros(npx,npy,128,8));
            phmat(:,:,:,1) = phasematsv(:,:,1:128);
            phmat(:,:,:,2) = phasematsv(:,:,123:250);
        end
    else
        if(pR.PhCo == 3)
            [~,~,nel] = size(phasematsv);
            AmpPhFull(:,:,jXducerSection+1) = pR.AmpPh1;
            phmat(:,:,1:nel,1+jXducerSection) = phasematsv;
        end
        pressure = pressure + pout;
        if test == 1
            pressure2 = pressure2 + pout2;
        end
    end
    pressar(:,:,:,jXducerSection) = pout;
    if(pR.PhCo == 3)
        fname=['AmpPh2_Z',num2str(jXducerSection),'_s',num2str(pR.json),'_c',num2str(pHAS.patnum),'.mat'];
        AmpPh1 = pR.AmpPh1;
        %         save(fname,'AmpPh1');
    end
end
if(pR.PhCo == 3)
    AmpPhar = [squeeze(AmpPhFull(:,:,1));squeeze(AmpPhFull(:,:,6));squeeze(AmpPhFull(:,:,2));squeeze(AmpPhFull(:,:,3)); ...
        squeeze(AmpPhFull(:,:,7));squeeze(AmpPhFull(:,:,8));squeeze(AmpPhFull(:,:,4));squeeze(AmpPhFull(:,:,5))];
    %fname = (['AmpPhar-PatNum',num2str(pHAS.patnum),'son',num2str(pR.json),'vel',num2str(pHAS.c(7)),'.mat']);
    %save(fname,'AmpPhar'); zach j dont save
else
    AmpPhar = [];
end
%    zpa1 = [ampphar(1:128,:);ampphar(257:384,:)];
%     zpa2 = ampphar(385:512,:);
%     zpa3 = ampphar(769:896,:);
%     zpa4 = ampphar(897:1024,:);
%     zpa5 = ampphar(129:256,:);
%     zpa6 = ampphar(513:640,:);
%     zpa7 = ampphar(641:768,:);
if test ==1
    % plotPhCo; 
    % Q_rel2=abs((pressure2.*conj(pressure2)).*aabs(Modl_in)*1e2*pERFA.fMHz./Z_modl);
    % cd([pR.pathOut,'/',pR.pname]);
    % save('Q_rel_noPhCo_fist_sonication.mat','Q_rel2');
end
%Q_rel=abs((pressure.*conj(pressure)).*absmodl_modl./Z_modl);        %Note the absmodl should be pressure attenuation DLP 10/12/20
Modl_Q = prepHAS_Model(Modl_in,0, 0,Foc_ijk,-pHAS.xtilt);
Q_rel=abs((pressure.*conj(pressure)).*aabs(Modl_in)*1e2*pERFA.fMHz./Z_modl);%%DLP 1/17/22 no need for absmodl_modl for just one use. save some memory %
%Q_rel=abs((pressure.*absmodl_modl.*conj(pressure.*absmodl_modl))./Z_modl);        %Note the absmodl should be pressure attenuation DLP 10/12/20
%q_apl=(pressure.*conj(pressure)).*1*1e2*pERFA.fMHz./Z_modl;
end

