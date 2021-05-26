% Visualization of arbitrary spherical harmonics transformations
%
% For details, see paper.
%
% Input:
%   T           SH-domain operation matrix of dimension (N_out+1)^2 x (N_in+1)^2
%   inputFormat (optional, default: 'N3D') Can be either 'N3D' or 'SN3D'
%               dependent on the spherical harmonics normalization scheme.
%               Default is N3D, but some audio folks prefer a different
%               normalization (--> ambix is SN3D)
%
%
%   Example 1: Rotation
%       load('examples/example_rotation.mat','T');
%       visualizeTransformation(T);
%
%   Example 2: Noise Reduction (direction-preserving multi-channel Wiener filter)
%       load('examples/example_noise_reduction_dp.mat','T');
%       visualizeTransformation(T);
%
%   Example 3: Noise Reduction (classic multi-channel Wiener filter)
%       load('examples/example_noise_reduction_pm.mat','T');
%       visualizeTransformation(T);
%
%   Example 4: Sound Field Translation using Adaptive Space Warping
%       load('examples/example_adaptive_space_warping.mat','T');
%       visualizeTransformation(T);
%
%
%--------------------------------------------------------------------------
% (c) 2021 - RWTH Aachen University
%--------------------------------------------------------------------------
% Version history:
% 1.0  - initial version - Maximilian Kentgens (kentgens@iks.rwth-aachen.de)
% 1.0a - standalone version which does not require SASP framework - Maximilian Kentgens (kentgens@iks.rwth-aachen.de)
%--------------------------------------------------------------------------
function hAx = visualizeTransformation(T,inputFormat)

    %% Check 
    if nargin<1
        error('Must provide operation matrix.');
    end
    
    %% Default spherical harmonics normalization scheme is N3D. 
    if nargin>=2 && ~isempty(inputFormat)
        switch lower(inputFormat)
            case 'n3d'
                T = T;
            case 'sn3d'
                inputOrders = floor(sqrt(0:size(T,2)-1));
                outputOrders = floor(sqrt(0:size(T,1)-1));
                T = diag(sqrt(2*outputOrders+1)) * T * diag(1./sqrt(2*inputOrders+1));
            otherwise
                error('Input format can be either N3D or SN3D.');
        end
    end
    
    %% Create Figure
    hFig = figure('Position',1*[0,0,1200,600]);
    hAx = axes;
    
    %% Draw everything
    visualizeDirectionalGain(T,hAx);
    visualizeEnergyVectors(T,hAx);
    
    %% Final axes formatting
    xlim([-pi,pi]);
    ylim([0,pi]);
    
    xticks([-pi -pi/2 0 pi/2 pi]);
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    yticks([0 pi/2 pi]);
    yticklabels({'0','\pi/2','\pi'});
    grid on;
    
    createAdditionalColorbar(hFig);
end

function visualizeDirectionalGain(T,hAx)
    [magnitude,magnitudeGrid] = calculateMagnitudes(T);
    p = DirectivityPlot();
    p.setAxes(hAx)
    p.setScatteredData(magnitudeGrid,magnitude,[],'linear');
    p.draw();
    hold on;
end

function visualizeEnergyVectors(T,hAx)
    
    directionsBeforeProcessing = getSphericalSamplingGrid('low');
    [directionsAfterProcessing, rE_coarse] = calculateDirectionsAfterProcessing(T,directionsBeforeProcessing);
    
    thetaBeforeProcessing = getCartGridInclination(directionsBeforeProcessing);
    phiBeforeProcessing   = getCartGridAzimuth(directionsBeforeProcessing);
    
    thetaAfterProcessing = getCartGridInclination(directionsAfterProcessing);
    phiAfterProcessing   = getCartGridAzimuth(directionsAfterProcessing);
    
    for dirIdx = 1:size(directionsBeforeProcessing,1)
        plotSlerp(directionsBeforeProcessing(dirIdx,:), ...
            directionsAfterProcessing(dirIdx,:), ...
            valueToColor(rE_coarse(dirIdx),[0,1]));
    end
    scatter(phiBeforeProcessing,thetaBeforeProcessing,'b','filled')
    scatter(phiAfterProcessing,thetaAfterProcessing,'r','filled')
end


function phi = getCartGridAzimuth(cartGrid)
    phi = atan2(cartGrid(:,2),cartGrid(:,1));
end

function theta = getCartGridInclination(cartGrid)
    r = sqrt(sum(cartGrid.^2,2));
    theta = acos(cartGrid(:,3)./r);
    theta(r==0) = 0;
end

function dist = getAngularDistance(cartGrid1,cartGrid2)
    r1 = sqrt(sum(cartGrid1.^2,2));
    r2 = sqrt(sum(cartGrid2.^2,2));
    tmp = sum(cartGrid1 .* cartGrid2,2) / (r1 .* r2);
    dist = acos(tmp);
    dist = real(dist); % DON'T KNOW WHY THIS IS NECESSARY...
end


function createAdditionalColorbar(hFig)
    hExistingCbr = colorbar();
    
    
    % dummy plot to create another colorbar
    hDummyFig = figure;
    hDummyAx = axes;
    hDummyPlot = plot([-pi,pi],[0,pi]);
    set(hDummyPlot,'Visible','off');
    axis(hDummyAx,'equal','tight');
    axis(hDummyAx,'off');
    hTmpCb = colorbar('AxisLocation','in');
    copiedObject = copyobj([hDummyAx,hTmpCb],hFig);
    hAdditiobalCbr = copiedObject(2);
    delete(hDummyFig);
    
    
    cbrPos = get(hExistingCbr,'Position');
    set(hExistingCbr,'Position',cbrPos+[0.07,0,-0.005,0])
    set(hAdditiobalCbr,'Position',cbrPos+[0.07,0,-0.005,0])
end


    
function [magnitudeRatio,sphericalDesign] = calculateMagnitudes(T)
    Nin = sqrt(size(T,2))-1;
    
    sphericalDesign = getSphericalSamplingGrid('medium');
    panningMatrix   = createSphericalHarmonicsMatrixCart(sphericalDesign,Nin)';
    
    magnitudeRatio = zeros(size(sphericalDesign,1),1);
    magnitudeBefore = sqrt(sum(panningMatrix(:,1).^2));
    for dirIdx = 1:size(sphericalDesign,1)
        coeffsAfter = T*panningMatrix(:,dirIdx);
        magnitudeAfter = sqrt(sum(coeffsAfter.^2));
        magnitudeRatio(dirIdx) = magnitudeAfter ./ magnitudeBefore;
    end
end


function plotSlerp(dir1,dir2,color)
    stepsPerDegree = 1.0;
    nTrajectory = 1 + round(stepsPerDegree * getAngularDistance(dir1,dir2) * 180/pi);
    trajectoryParameter = linspace(0,1,nTrajectory);
    cartTrajectory = zeros(nTrajectory,3);
    for trajectoryIdx = 1:nTrajectory
        alpha = trajectoryParameter(trajectoryIdx);
        cartTrajectory(trajectoryIdx,:) = alpha*dir1 + (1-alpha)*dir2;
    end
    
    thetaTrajectory = getCartGridInclination(cartTrajectory);
    phiTrajectory = getCartGridAzimuth(cartTrajectory);
    
    [thetaTrajectory,phiTrajectory] = repairTrajectoryWrapAround(thetaTrajectory,phiTrajectory);
    
    line(phiTrajectory,thetaTrajectory,'Color',[0,0,0],'LineWidth',3);
    line(phiTrajectory,thetaTrajectory,'Color',color,'LineWidth',2);
end

function [thetaTrajectory,phiTrajectory] = repairTrajectoryWrapAround(thetaTrajectory,phiTrajectory)
    diffPhi = diff(phiTrajectory);
    a = abs(diffPhi) > abs(diffPhi+pi);
    b = abs(diffPhi) > abs(diffPhi-pi);
    if any(a) || any(b)
        idx = union(find(a,1),find(b,1));
        sign = 1.0*any(a) - 1.0*any(b);
        phiInsert = [phiTrajectory(idx+1)+2*pi*sign;NaN;phiTrajectory(idx)-2*pi*sign];
        thetaInsert = [thetaTrajectory(idx+1);NaN;thetaTrajectory(idx)];
        phiTrajectory = [phiTrajectory(1:idx); phiInsert; phiTrajectory(idx+1:end)];
        thetaTrajectory = [thetaTrajectory(1:idx); thetaInsert; thetaTrajectory(idx+1:end)];
    end
end

function [directionsAfterProcessing, rE] = calculateDirectionsAfterProcessing(T,directionsBeforeProcessing)
    Nin = sqrt(size(T,2))-1;
    panningMatrix = createSphericalHarmonicsMatrixCart(directionsBeforeProcessing,Nin)';
    directionsAfterProcessing = zeros(size(directionsBeforeProcessing));
    rE = zeros(size(directionsBeforeProcessing,1),1);
    for dirIdx = 1:size(directionsBeforeProcessing,1)
        inputCoeffs = panningMatrix(:,dirIdx);
        outputCoeffs = T * inputCoeffs;
        energyVector = calculateEnergyVector(outputCoeffs);
        rE(dirIdx) = sqrt(energyVector(:)'*energyVector(:));
        directionsAfterProcessing(dirIdx,:) = energyVector(:)' ./ norm(energyVector);
    end
end


function rgb = valueToColor(values, colorRange)
    cmap = colormap();
    nColors = size(cmap,1);
    %assert(size(cmap,1)==256)
    % Normalize the values to be between 1 and 256
    values(values < colorRange(1)) = colorRange(1);
    values(values > colorRange(2)) = colorRange(2);
    valsN = round(((values - colorRange(1)) ./ diff(colorRange)) .* (nColors-1))+1;
    % Convert any nans to ones
    valsN(isnan(valsN)) = 1;
    % Convert the normalized values to the RGB values of the colormap
    rgb = cmap(valsN, :);
end


% Calculate energy vector for given SH coefficiens
function rE = calculateEnergyVector(coeffs)
    N = sqrt(numel(coeffs))-1;

    % Only calculate spherical harmonics matrix once and store in memory
    persistent storeN storeCartSphGrid storeBeamMatrix
    if isempty(storeN) || storeN~=N
        storeN            = N;
        storeCartSphGrid  = getSphericalSamplingGrid('high');
        storeBeamMatrix   = createSphericalHarmonicsMatrixCart(storeCartSphGrid,storeN);
    end
    
    % Approximate integral over unit sphere for calculation of energy
    % vector by summing over a dense spherical grid
    sigSpatial = storeBeamMatrix * coeffs(:);
    denominator = sum(abs(sigSpatial).^2);
    rE_x = sum(storeCartSphGrid(:,1) .* abs(sigSpatial).^2) / denominator;
    rE_y = sum(storeCartSphGrid(:,2) .* abs(sigSpatial).^2) / denominator;
    rE_z = sum(storeCartSphGrid(:,3) .* abs(sigSpatial).^2) / denominator;
    rE = [rE_x, rE_y, rE_z];
    
    % Basic sanity check
    assert(norm(rE)<=1);
end

% Create spherical harmonics matrix with Cartesian grid input
function Y = createSphericalHarmonicsMatrixCart(cartGrid,shOrder)
    thetaVec = getCartGridInclination(cartGrid);
    phiVec = getCartGridAzimuth(cartGrid);
    Y = createSphericalHarmonicsMatrix(thetaVec,phiVec,shOrder);
end



% Return spherical grid design of requested resolution.
%
% Return value is matrix Qx3 where Q is the number of sampling points. Each
% row of the matrix is a grid point in Cartesian coordinates (x,y,z).
function cartGrid = getSphericalSamplingGrid(type)
    persistent cartGridHighDensity
    persistent cartGridMediumDensity
    persistent cartGridLowDensity
        
    %%
    % load grids from file 'spherical_sampling_grids.mat' if
    % not already in memory
    if isempty(cartGridHighDensity) ...
            || isempty(cartGridMediumDensity) ...
            || isempty(cartGridLowDensity)
        [directory,~,~] = fileparts(mfilename('fullpath'));
        load(fullfile(directory,'spherical_sampling_grids.mat'),...
            'cartGridHighDensity', ...
            'cartGridMediumDensity', ...
            'cartGridLowDensity');
        % perform some basic data sanity checks
        assert(size(cartGridHighDensity,2)==3);
        assert(size(cartGridMediumDensity,2)==3);
        assert(size(cartGridLowDensity,2)==3);
        assert(size(cartGridLowDensity,1)<size(cartGridMediumDensity,1));
        assert(size(cartGridMediumDensity,1)<size(cartGridHighDensity,1));
    end
    
    %% return requested grid
    switch lower(type)
        case 'high'
            cartGrid = cartGridHighDensity;
        case 'medium'
            cartGrid = cartGridMediumDensity;
        case 'low'
            cartGrid = cartGridLowDensity;
        otherwise
            assert(false);
    end
end