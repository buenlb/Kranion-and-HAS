function [Vhas, pHAS] = addCtModel(boneIndx,Vhas,CT,pHAS)

USE_RAW_HU = 1;

HUbone = 3500;

if ~USE_RAW_HU
    %% Directly from Leung, assuming energy is matched
    f = CT/HUbone;
    porosity = 100*(1-f);
    
    binWidth = 10;
    bins = 0:binWidth:100+binWidth/2;
    binnedF = -(bins-100)/100;
    binnedHu = binnedF*HUbone;
    atten = zeros(size(bins));
    
    atten(bins<=25) = (0.522-0.922)/(0-25)*bins(bins<=25)+0.522;
    atten(bins<=40 & bins>25) = (0.922-2.136)/(25-40)*bins(bins<=40 & bins>25)+0.922-(0.922-2.136)/(25-40)*25;
    atten(bins<=60 & bins>40) = 2.136;
    atten(bins<=85 & bins>60) = (2.136-0.1183)/(60-85)*bins(bins<=85 & bins>60)+2.136-(2.136-0.1183)/(60-85)*60;
    atten(bins<=100 & bins>85) = (0.1183-0.0412)/(85-100)*bins(bins<=100 & bins>85)+0.1183-(0.1183-0.0412)/(85-100)*85;
    
    rho = 1920*binnedF+1e3*(1-binnedF);
    c = 0.75*binnedHu+1320;
    
    binnedCT = zeros(size(CT));
    for ii = 1:length(bins)
        binnedCT(porosity >= bins(ii)-binWidth/2 & porosity < bins(ii)+binWidth/2) = ii;
    end
    binnedCT(porosity==100) = ii;
    
    for ii = 1:length(boneIndx)
        Vhas(Vhas==boneIndx(ii)) = min(boneIndx);
    end
    
    pHAS.a = [pHAS.a(1:(min(boneIndx)-1)),atten]/0.65; % The 0.65 converts this from Np/cm to Np/cm/MHz
    pHAS.a_abs = [pHAS.a_abs(1:(min(boneIndx)-1)),ones(size(atten))*0.47];
    pHAS.rho = [pHAS.rho(1:(min(boneIndx)-1)),rho];
    pHAS.c = [pHAS.c(1:(min(boneIndx)-1)),c];
    pHAS.randvc = zeros(size(pHAS.c));
    
    Vhas(Vhas == min(boneIndx)) = binnedCT(Vhas == min(boneIndx))+min(boneIndx)-1;
else
    huVerts = [0,300,800,1200,1500,2000];
    atten_verts = [0.041, 0.118, 2.136, 2.136, 0.922, 0.522];
    hu = 0:2000;

    % Interpolate attenuation values for the given HU range
    atten = interp1(huVerts, atten_verts, hu, 'linear', 'extrap');
    rho = 0.46*hu+1e3;
    c = 0.75*hu+1320;

    for ii = 1:length(boneIndx)
        Vhas(Vhas==boneIndx(ii)) = min(boneIndx);
    end

    pHAS.a = [pHAS.a(1:(min(boneIndx)-1)),atten]/0.65; % The 0.65 converts this from Np/cm to Np/cm/MHz
    pHAS.a_abs = [pHAS.a_abs(1:(min(boneIndx)-1)),ones(size(atten))*0.47];
    pHAS.rho = [pHAS.rho(1:(min(boneIndx)-1)),rho];
    pHAS.c = [pHAS.c(1:(min(boneIndx)-1)),c];
    pHAS.randvc = zeros(size(pHAS.c));

    binnedCT = round(CT);
    binnedCT(binnedCT<0) = 0;
    binnedCT(binnedCT>max(hu))=max(hu);
    Vhas(Vhas == min(boneIndx)) = binnedCT(Vhas == min(boneIndx))+min(boneIndx);
end