% Visualization class for directivity plots
%
%   Example 1:
%       p = DirectivityPlot();
%       p.plotType = 'spherical';
%       %p.plotType = 'lobes';
%       %p.plotType = 'equirectangular';
%       axisPhi = linspace(-pi,pi,144);
%       axisTheta = linspace(0,pi,72);
%       data = sin(axisTheta(:)).^2*cos(2*axisPhi(:)');
%       p.setData(axisTheta,axisPhi,data);
%       p.draw();
%
%   Example 2: Scattered data interpolation
%       coordx = [ ...
%           0.8506         0   -0.5258;
%           0.5257   -0.8506   -0.0000;
%           -0.0000   -0.5258    0.8506;
%           0.8506         0    0.5257;
%           -0.5258   -0.8506   -0.0000;
%           -0.0000    0.5257   -0.8506;
%           -0.8506   -0.0000   -0.5258;
%           -0.5258    0.8506   -0.0000;
%           -0.0000    0.5258    0.8506;
%           -0.8506   -0.0000    0.5257;
%           0.5257    0.8506   -0.0000;
%           -0.0000   -0.5257   -0.8506 ];
%       scattered_data = rand(size(coordx,1),1);
%       p = DirectivityPlot();
%       p.setScatteredData(coordx,scattered_data,[],'linear');
%       p.draw();
%
%--------------------------------------------------------------------------
% (c) 2021 - RWTH Aachen University
%--------------------------------------------------------------------------
% Version history:
% 1.0  - initial version - Maximilian Kentgens (kentgens@iks.rwth-aachen.de)
% 1.0a - standalone version which does not require SASP framework - Maximilian Kentgens (kentgens@iks.rwth-aachen.de)
%--------------------------------------------------------------------------
classdef DirectivityPlot < handle

    properties
        plotType = 'equirectangular';
        data
        axisPhi
        axisTheta
        bgPicture
        bgAlpha
        display = 'linear'
    end
    
    properties(Access=private)
        hPlot
        hBackground
        hAxes = []
    end
    
    methods
        function this = DirectivityPlot()
        end
        
        function setData(this,axisTheta,axisPhi,data)
            axisPhi(axisPhi>pi) = axisPhi(axisPhi>pi) - 2*pi;
            this.axisTheta = axisTheta;
            this.axisPhi   = axisPhi;
            this.data      = data;
        end
        
        function setScatteredData(this,coordx,scattered_data,resolution,method)
            if nargin<4 || isempty(resolution)
                resolution = pi/180;
            end
            if nargin<5 || isempty(method)
                method = 'linear'; % linear, nearest
            end
            %% create grid
            axisPhi = -pi:resolution:pi;
            axisTheta = 0:resolution:pi;
            [phix,thetax] = meshgrid(axisPhi,axisTheta);
            %% grid in cartesian coordinates
            xx = sin(thetax) .* cos(phix);
            yy = sin(thetax) .* sin(phix);
            zz = cos(thetax);
            %% do the interpolation in cartesian coordinates
            if isa(coordx,'SASP.Coordinate')
                grid = [[coordx.x]', [coordx.y]', [coordx.z]'];
            else
                % asumme data to be cartesian coordinate matrix of
                % dimension Qx3
                grid = coordx;
            end
            grid = grid ./ sqrt(sum(grid.^2,2));
            F = scatteredInterpolant(grid(:,1),grid(:,2),grid(:,3),scattered_data(:),method);
            %% convert the data to logarithmic represantation
            if strcmp(this.display, 'linear')
                data = F(xx,yy,zz);
            else
                data = 20*log(abs(F(xx,yy,zz)));
            end
            %% save data
            this.setData(axisTheta,axisPhi,data);
        end
        
        function setBackgroundPicture(this,filename,color,alpha,hist_limits_from,hist_limit_to)
            if ~isempty(filename)
                I = imread(filename);
                if nargin>=3 && ~color
                    I = rgb2gray(I);
                    I = cat(3,I,I,I);
                end
                if nargin<4
                    alpha = 0.5;
                end
                if nargin>=6 && ~isempty(hist_limits_from) && ~isempty(hist_limits_to)
                    I = imadjust(I,hist_limits_from,hist_limit_to);
                end
            else
                I = [];
                alpha = [];
            end
            this.bgPicture = I;
            this.bgAlpha   = alpha;
        end
        
        function draw(this)
            if isempty(this.hAxes) || ~ishghandle(this.hAxes)
                figure;
                this.hAxes = axes;
            end
            switch lower(this.plotType)
                case 'equirectangular'
                    this.drawEquirectangular();
                case 'lobes'
                    this.drawLobes();
                case 'spherical'
                    this.drawSpherical();
                otherwise
                    error('Unknown plot type specified.');
            end
        end
        
        function setAxes(this,ax)
            this.hAxes = ax;
        end
    end
    
    methods(Access=private)
        function drawEquirectangular(this)
            if isempty(this.hPlot) || ~ishghandle(this.hPlot) || ~strcmpi(this.hPlot.UserData,this.plotType)
                cla(this.hAxes);
                if ~isempty(this.bgPicture)
                    xx = linspace(-pi,pi,size(this.bgPicture,2));
                    yy = linspace(0,pi,size(this.bgPicture,1));
                    this.hBackground = imagesc(this.hAxes,xx,yy,fliplr(this.bgPicture));
                    hold(this.hAxes,'on');
                end
                this.hPlot = imagesc(this.hAxes,this.axisPhi,this.axisTheta,this.data);
                this.hPlot.UserData = this.plotType;
                if ~isempty(this.bgPicture)
                    alpha(this.hPlot,this.bgAlpha);
                end
                this.decorateEquirectangular();
            else
                this.hPlot.XData    = this.axisPhi;
                this.hPlot.YData    = this.axisTheta;
                this.hPlot.CData    = this.data;
            end
        end
        
        function decorateEquirectangular(this)
            colorbar(this.hAxes);
            %axis(this.hAxes,'ij');
            axis(this.hAxes,'equal');
            axis(this.hAxes,'tight');
            set(this.hAxes,'YDir','reverse');
            set(this.hAxes,'XDir','reverse');
            ylabel(this.hAxes,'inclination \theta [rad]');
            xlabel(this.hAxes,'azimuth \phi [rad]');
        end
        
        function drawSpherical(this)
            [xx,yy,zz] = this.sphereEx(this.axisTheta,this.axisPhi);
            if isempty(this.hPlot) || ~ishghandle(this.hPlot) || ~strcmpi(this.hPlot.UserData,this.plotType)
                cla(this.hAxes);
                [warpx,warpy,warpz] = sphere(30);
                if ~isempty(this.bgPicture)
                    this.hBackground = warp(warpx*0.99,warpy*0.99,warpz*0.99,flipud(this.bgPicture));
                    hold(this.hAxes,'on');
                end
                this.hPlot = surf(this.hAxes,xx,yy,zz,this.data,'EdgeColor','none');
                this.hPlot.UserData = this.plotType;
                if ~isempty(this.bgPicture)
                    alpha(this.hPlot,this.bgAlpha);
                end
                this.decorateSpherical();
            else
                this.hPlot.XData    = xx;
                this.hPlot.YData    = yy;
                this.hPlot.ZData    = zz;
                this.hPlot.CData    = this.data;
            end
        end
        
        function [xx,yy,zz] = sphereEx(this,theta,phi)
            theta = theta(:);
            phi   = phi(:)';
            xx = sin(theta) * cos(phi);
            yy = sin(theta) * sin(phi);
            zz = cos(theta) * ones(size(phi));
        end
        
        function decorateSpherical(this)
            axis(this.hAxes,'equal');
            colorbar(this.hAxes);
            xlabel(this.hAxes,'x');
            ylabel(this.hAxes,'y');
            zlabel(this.hAxes,'z');
        end
        
        function drawLobes(this)
            [xx,yy,zz] = this.sphereEx(this.axisTheta,this.axisPhi);
            xx = xx .* abs(this.data);
            yy = yy .* abs(this.data);
            zz = zz .* abs(this.data);
            if isempty(this.hPlot) || ~ishghandle(this.hPlot) || ~strcmpi(this.hPlot.UserData,this.plotType)
                %cla(this.hAxes);
                %this.hPlot = surf(this.hAxes,xx,yy,zz,sign(this.data));
                this.hPlot = surf(this.hAxes,xx,yy,zz,this.data, 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
                this.hPlot.UserData = this.plotType;
                this.decorateLobes();
            else
                this.hPlot.XData    = xx;
                this.hPlot.YData    = yy;
                this.hPlot.ZData    = zz;
                this.hPlot.CData    = this.data;
            end
        end
        
        function decorateLobes(this)
            axis(this.hAxes,'equal');
            xlabel(this.hAxes,'x');
            ylabel(this.hAxes,'y');
            zlabel(this.hAxes,'z');
        end
        
    end
end

