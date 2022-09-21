function [apertures app_lim centrex centrey] = make_apertures(imin,spacing,radius,sigma,string,dither,xedge,yedge)
% Function to create randomly or ordered distribution of aperures for
% Ptychography
% imin = input image
% nApert =  desired number of apertures
% radius = radius of aperture
% sigma = sigma of gaussian
% string = string defining aperture type (grid, spiral or square)

%% Check the size of the image (2D or 3D)
if size(size(imin),2) > 2
    [imy imx imz] = size(imin);
else
    [imy imx] = size(imin);
end

xfar = (imx - xedge)+2;
yfar = (imy - yedge)+2;
%% begin to construct Gaussian function at each appeture
half_size = ceil(imy/2);

diameter = 2*radius;

if strcmp(string,'spiral') == 1
    
    [rr cc] = meshgrid(1:imy); % coordinates to define where apertures are
    
    
    r_lin = reshape(rr, 1,imy^2);
    c_lin = reshape(cc, 1,imy^2);
    r_lin(r_lin < xedge) = -1;
    r_lin(r_lin > xfar) = -1;
    c_lin(c_lin < yedge) = -1;
    c_lin(c_lin > yfar) = -1;
%     centre_list = [r_lin;c_lin];
%     count = 1;
    
        xlim = size(imin,1); % size in 1st dimension
    ylim = size(imin,2); % size in 2nd dimension
    
    xcent = ceil(xlim/2); % centre of image in first dimension
    ycent = ceil(ylim/2); % centre of the image in 2nd dimension
    
%     [x,y]=meshgrid(1:xlim,1:ylim);
    x = reshape(r_lin,[imy, imy]);
    y = reshape(c_lin,[imy,imy]);
    [angle,r]=cart2pol(x-xcent,y-ycent); % Radius and angle used to define spiral
    
    Spiral=sin(r/(12*pi)); % + angle); %/(4*pi) /(2*pi)
    %     imagesc(Spiral);
    %     axis equal; axis off;
    
    ind = find(Spiral > 0.9); % 0.9 Find the peaks of the Spiral
    
    places = 1:1:size(ind,1);
    
    [centrex_in, centrey_in] = ind2sub([xlim,ylim],ind);
    
    [rr cc] = meshgrid(1:ylim);
    count = 1;
    for i = 1:size(centrex_in,1)
        if mod(centrex_in(i),spacing) == 0 && mod(centrey_in(i),spacing) == 0 && x(centrex_in(i),centrey_in(i)) ~= -1 && y(centrex_in(i),centrey_in(i)) ~= -1
            theta_r = pi/3;
            tmp_x = 3*cos(theta_r)*(rr - centrex_in(i)) - 3*sin(theta_r)*(cc - centrey_in(i));
            tmp_y = sin(theta_r)*(rr - centrex_in(i)) + cos(theta_r)*(cc - centrey_in(i));
%             tmp_xy = [rr - centrex_in(i) cc - centrey_in(i)]*[3*cos(theta_r) -sin(theta_r); 3*sin(theta_r) cos(theta_r)];
            modulus = sqrt((tmp_x).^2+(tmp_y).^2);
            %             modulus = sqrt((rr - centrex_in(i)).^2+(cc - centrey_in(i)).^2);
            gauss_app = exp(- modulus.^2/(2*sigma^2)); % create 2D gaussian
            gauss_app=gauss_app/max(max(gauss_app)); % normalise to 1
            centrex(count) = centrex_in(i);
            centrey(count) = centrey_in(i);
            if size(size(imin),2) > 2
                app_lim(:,:,:,i/spacing) = repmat(sqrt((rr - centrex_in(i)).^2+(cc - centrey_in(i)).^2)<=radius,[1 1 imz]);
                apertures(:,:,:,i/spacing) = app_lim(:,:,:,i/spacing).*repmat(gauss_app,[1 1 imz]); % create truncated gaussian aperture
            else
                app_lim(:,:,count) = sqrt((rr - centrex_in(i)).^2+(cc - centrey_in(i)).^2)<=radius;
                apertures(:,:,count) = app_lim(:,:,count).*gauss_app; %fresnel_advance(app_lim(:,:,i/spacing),1,1,5,1); % create truncated gaussian aperture
            end
            count = count + 1;
        end
        
    end
    
    
%     xlim = size(imin,1); % size in 1st dimension
%     ylim = size(imin,2); % size in 2nd dimension
%     
%     xcent = ceil(xlim/2); % centre of image in first dimension
%     ycent = ceil(ylim/2); % centre of the image in 2nd dimension
%     
%     [x,y]=meshgrid(1:xlim,1:ylim);
%     [angle,r]=cart2pol(x-xcent,y-ycent); % Radius and angle used to define spiral
%     
%     Spiral=sin(r + angle);
%     %     imagesc(Spiral);
%     %     axis equal; axis off;
%     
%     ind = find(Spiral > 0.9); % Find the peaks of the Spiral
%     
%     places = 1:4:size(ind,1);
%     
%     [centrex, centrey] = ind2sub([xlim,ylim],ind);
%     
%     [rr cc] = meshgrid(1:ylim);
%     
%     for i = 1:size(centrex,1)
%         if mod(i,spacing) == 0
%             modulus = sqrt((rr - centrex(i)).^2+(cc - centrey(i)).^2);
%             gauss_app = exp(- modulus.^2/(2*sigma^2)); % create 2D gaussian
%             gauss_app=gauss_app/max(max(gauss_app)); % normalise to 1
%             if size(size(imin),2) > 2
%                 app_lim(:,:,:,i/spacing) = repmat(sqrt((rr - centrex(i)).^2+(cc - centrey(i)).^2)<=radius,[1 1 imz]);
%                 apertures(:,:,:,i/spacing) = app_lim(:,:,:,i/spacing).*repmat(gauss_app,[1 1 imz]); % create truncated gaussian aperture
%             else
%                 app_lim(:,:,i/spacing) = sqrt((rr - centrex(i)).^2+(cc - centrey(i)).^2)<=radius;
%                 apertures(:,:,i/spacing) = app_lim(:,:,i/spacing).*gauss_app; %fresnel_advance(app_lim(:,:,i/spacing),1,1,5,1); % create truncated gaussian aperture
%             end
%         end
%     end
    
elseif strcmp(string, 'grid') == 1
    [rr cc] = meshgrid(1:imy); % coordinates to define where apertures are
    
    
    r_lin = reshape(rr, 1,imy^2);
    c_lin = reshape(cc, 1,imy^2);
    r_lin(r_lin < xedge) = -1;
    r_lin(r_lin > xfar) = -1;
    c_lin(c_lin < yedge) = -1;
    c_lin(c_lin > yfar) = -1;
    centre_list = [r_lin;c_lin];
    count = 1;
    for i = 1:size(centre_list,2)
        if mod(centre_list(1,i),spacing) == 0 && mod(centre_list(2,i),spacing) == 0 && centre_list(1,i) ~= -1 && centre_list(2,i) ~= -1
            ditherx = round(rand(1)*dither - (dither/2));
            dithery = round(rand(1)*dither - (dither/2));
            centrex(count) = (centre_list(1,i)+ditherx);
            centrey(count) = (centre_list(2,i) + dithery);
            
%             modulus = sqrt((rr - centrex(count)).^2+(cc - centrey(count)).^2);
%             gauss_app = exp(- modulus.^2/(2*sigma^2)); % cr   ditherx = round(rand(1)*dither - (dither/2));
%             dithery = round(rand(1)*dither - (dither/2));
%             centrex(count) = (centre_list(1,i)+ditherx);
%             centrey(count) = (centre_list(2,i) + dithery);
%             
%             modulus = sqrt((rr - centrex(count)).^2+(cc - centrey(count)).^2);

            
            theta_r = pi/3;
            tmp_x = 3*cos(theta_r)*(rr - centrex(count)) - 3*sin(theta_r)*(cc - centrey(count));
            tmp_y = sin(theta_r)*(rr - centrex(count)) + cos(theta_r)*(cc - centrey(count));
%             tmp_xy = [rr - centrex_in(i) cc - centrey_in(i)]*[3*cos(theta_r) -sin(theta_r); 3*sin(theta_r) cos(theta_r)];
            modulus = sqrt((tmp_x).^2+(tmp_y).^2);

            gauss_app = exp(- modulus.^2/(2*sigma^2)); % create 2D gaussianeate 2D gaussian
            gauss_app=gauss_app/max(max(gauss_app)); % normalise to 1
            if size(size(imin),2) > 2
                app_lim(:,:,:,count) = repmat(sqrt((rr - centrex(count)).^2+(cc - centrey(count)).^2)<=radius,[1 1 imz]);% define extent of apeture
                apertures(:,:,:,count) = app_lim(:,:,:,count).*repmat(gauss_app,[1 1 imz]);%repmat(sqrt((rr - centrex(count)).^2+(cc - centrey(count)).^2)<=radius,[1 1 imz]); % create truncated gaussian aperture
              
            else
                app_lim(:,:,count) = sqrt((rr - centrex(count)).^2+(cc - centrey(count)).^2)<=radius;
                apertures(:,:,count) = app_lim(:,:,count).*gauss_app;
            end
            count = count+1;
        end
    end
    
    
    
elseif strcmp(string, 'square') == 1
    [rr cc] = meshgrid(1:imy); % coordinates to define where apertures are
    
    
    r_lin = reshape(rr, 1,imy^2);
    c_lin = reshape(cc, 1,imy^2);
    r_lin(r_lin < xedge) = -1;
    r_lin(r_lin > xfar) = -1;
    c_lin(c_lin < yedge) = -1;
    c_lin(c_lin > yfar) = -1;
    centre_list = [r_lin;c_lin];
    count = 1;
    for i = 1:size(centre_list,2)
        if mod(centre_list(1,i),spacing) == 0 && mod(centre_list(2,i),spacing) == 0 && centre_list(1,i) ~= -1 && centre_list(2,i) ~= -1
            ditherx = round(rand(1)*dither - (dither/2));
            dithery = round(rand(1)*dither - (dither/2));
            centrex(count) = (centre_list(1,i)+ditherx);
            centrey(count) = (centre_list(2,i) + dithery);
            
            modulus = sqrt((rr - centrex(count)).^2+(cc - centrey(count)).^2);
            gauss_app = exp(- modulus.^2/(2*sigma^2)); % create 2D gaussian
            gauss_app=gauss_app/max(max(gauss_app)); % normalise to 1
            %                 if count == 1000
            %                 imagesc(gauss_app); axis image
            %                 centrex(count)
            %                 centrey(count)
            %                 pause()
            %                 end
            if size(size(imin),2) > 2
                app_holder = zeros(imx,imy);
                app_holder(centrey(count) - radius:centrey(count) + (radius - 1),...
                    centrex(count) - radius:centrex(count) + (radius - 1)) = 1;
                app_lim(:,:,:,count) = repmat(app_holder,[1 1 imz]); % define extent of apeture
                apertures(:,:,:,count) = app_lim(:,:,:,count).*repmat(gauss_app,[1 1 imz]); % create truncated gaussian aperture
                clear app_holder
            else
                app_holder = zeros(imx,imy);
                app_holder(centrey(count) - radius:centrey(count) + (radius - 1),...
                    centrex(count) - radius:centrex(count) + (radius - 1)) = 1;
                app_lim(:,:,count) = app_holder;
                apertures(:,:,count) = app_lim(:,:,count).*gauss_app;
                clear app_holder
            end
            count = count+1;
        end
    end
    
end


if strcmp(string, 'big') == 1
    [rr cc] = meshgrid(1:imy); % coordinates to define where apertures are
    
    
    r_lin = reshape(rr, 1,imy^2);
    c_lin = reshape(cc, 1,imy^2);
    r_lin(r_lin < xedge) = -1;
    r_lin(r_lin > xfar) = -1;
    c_lin(c_lin < yedge) = -1;
    c_lin(c_lin > yfar) = -1;
    centre_list = [r_lin;c_lin];
    count = 1;
    for i = 1:size(centre_list,2)
        if mod(centre_list(1,i),spacing) == 0 && mod(centre_list(2,i),spacing) == 0 && centre_list(1,i) ~= -1 && centre_list(2,i) ~= -1
            ditherx = round(rand(1)*dither - (dither/2));
            dithery = round(rand(1)*dither - (dither/2));
            centrex(count) = (centre_list(1,i)+ditherx);
            centrey(count) = (centre_list(2,i) + dithery);
            
            modulus = sqrt((rr - centrex(1)).^2+(cc - centrey(1)).^2);
            gauss_app = exp(- modulus.^2/(2*sigma^2)); % create 2D gaussian
            gauss_app=gauss_app/max(max(gauss_app)); % normalise to 1
            if size(size(imin),2) > 2
                app_lim(:,:,:) = repmat(sqrt((rr - centrex(1)).^2+(cc - centrey(1)).^2)<=radius,[1 1 imz]);% define extent of apeture
                apertures(:,:,:) = app_lim.*repmat(gauss_app,[1 1 imz]);%repmat(sqrt((rr - centrex(1)).^2+(cc - centrey(1)).^2)<=radius,[1 1 imz]); % create truncated gaussian aperture
              
            else
                app_lim(:,:,count) = sqrt((rr - centrex(count)).^2+(cc - centrey(count)).^2)<=radius;
                apertures(:,:,count) = app_lim(:,:,count).*gauss_app;
            end
            count = count+1;
        end
    end
end
end









    
    
    
    
    
    
    
    
    