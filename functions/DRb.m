%% sir-DR reconstruction algorithm

function [big_obj,aperture,fourier_error,initial_obj,initial_aperture, hist] = DRb(ePIE_inputs,varargin)
%varargin = {beta_obj, beta_ap, probeNorm, init_weight, final_weight, order, semi_implicit_P}
rng('shuffle','twister');
%% setup working and save directories
imout = 1; % boolean to determine whether to monitor progress or not

[~,jobID] = system('echo $JOB_ID');
jobID = jobID(~isspace(jobID));
if isempty(jobID)
    time_str = char(datetime('now','Format','uuuuMMdd_HHmm'));
    jobID = sprintf('local_%s',time_str);
end
%% Load inputs from struct
diffpats = ePIE_inputs(1).Patterns;
positions = ePIE_inputs(1).Positions;
pixel_size = ePIE_inputs(1).PixelSize;
big_obj = ePIE_inputs(1).InitialObj;
aperture_radius = ePIE_inputs(1).ApRadius;
aperture = ePIE_inputs(1).InitialAp;
iterations = ePIE_inputs(1).Iterations;

object = ePIE_inputs(1).object;
hist = struct('amp',[],'angle',[],'Rfactor',[]);

filename = ePIE_inputs(1).FileName;
filename = strcat('rPIE_',filename,'_',jobID);
filename = strrep(filename,'__','_');
probe_mask_flag = 0;

%% === input parameters === %%
optional_args = {.9, .1, 0.9, 0, 0.1, 0.4, 4, 0}; %default values for varargin parameters
nva = length(varargin);
optional_args(1:nva) = varargin;
[beta_obj,beta_ap,alpha, probeNorm, init_weight,final_weight,order, semi_implicit_P] = optional_args{:};

if isfield(ePIE_inputs, 'saveOutput'), saveOutput = ePIE_inputs.saveOutput;
else  saveOutput = 0;
end

if isfield(ePIE_inputs, 'save_path'), save_string = ePIE_inputs.save_path;
else save_string = [ pwd '/Results_ptychography/']; % Place to save results
end

if isfield(ePIE_inputs, 'updateAp'), update_aperture=ePIE_inputs.updateAp;
else update_aperture = 1;
end

if isfield(ePIE_inputs, 'GpuFlag'),  gpu = ePIE_inputs(1).GpuFlag;
else gpu = 0;
end

if isfield(ePIE_inputs, 'showim'),  showim = ePIE_inputs(1).showim;
else showim = 0;
end

if isfield(ePIE_inputs, 'do_posi'), do_posi = ePIE_inputs.do_posi;
else do_posi = 0;
end

if isfield(ePIE_inputs, 'ApType'), ApType = ePIE_inputs.ApType;
else ApType = 'R';
end


if isfield(ePIE_inputs, 'save_intermediate');
    save_intermediate = ePIE_inputs.save_intermediate;
else save_intermediate = 0;
end

if isfield(ePIE_inputs, 'save_inter_freq');
    save_inter_freq = ePIE_inputs.save_inter_freq;
else save_inter_freq = 100;
end

if isfield(ePIE_inputs, 'obj_scale'), obj_scale = ePIE_inputs.obj_scale;
else obj_scale = 1;
end

if isfield(ePIE_inputs, 'freeze_aperture');
    freeze_aperture = ePIE_inputs.freeze_aperture;
else freeze_aperture = Inf;
end

if isfield(ePIE_inputs, 'R_probe_mask'); %real space probe mask
    R_probe_mask = single(ePIE_inputs.R_probe_mask);
    %params.R_probe_mask = R_probe_mask;
    probe_mask_flag = 1;
end

if isfield(ePIE_inputs, 'F_probe_mask'); %k-space probe mask (don't fftshift)
    F_probe_mask = single(ePIE_inputs.F_probe_mask);
    %params.F_probe_mask = F_probe_mask;
    probe_mask_flag = 2;
end

%% print parameters
fprintf('filename = %s\n', filename);
fprintf('iterations = %d\n', iterations);
fprintf('gpu flag = %d\n', gpu);
fprintf('updating probe = %d\n', update_aperture);
fprintf('beta probe = %0.5f\n', beta_ap);
fprintf('beta object = %0.2f\n', beta_obj);
fprintf('probe norm = %d\n', probeNorm);
fprintf('Initial weight = %.2f, final weight = %.2f, order = %d\n',init_weight, final_weight, order)
fprintf('semi implicit on P = %d\n',semi_implicit_P);
fprintf('probe_mask_flag = %d\n',probe_mask_flag);
fprintf('Number of diff patterns = %d\n',size(diffpats,3));
fprintf('big_obj size = %dx%d, aperture size = %dx%d, pixel_size = %f\n',size(big_obj), size(aperture),pixel_size);
clear ePIE_inputs
%% Define parameters from data and for reconstruction
for ii = 1:size(diffpats,3)
    diffpats(:,:,ii) = fftshift(sqrt(diffpats(:,:,ii)));
end
diffpats = single(diffpats);
[y_kspace,~] = size(diffpats(:,:,1)); % Size of diffraction patterns
[N1,N2,nApert] = size(diffpats);
little_area = y_kspace; % Region of reconstruction for placing back into full size image

%% Get centre positions for cropping (should be a 2 by n vector)
%positions = convert_to_pixel_positions(positions,pixel_size);
[pixelPositions, bigx, bigy] = ...
    convert_to_pixel_positions_testing5(positions,pixel_size,little_area);
centrey = round(pixelPositions(:,2));
centrex = round(pixelPositions(:,1));
%centBig = round((bigx+1)/2);
Y1 = centrey - floor(y_kspace/2); Y2 = Y1+y_kspace-1;
X1 = centrex - floor(y_kspace/2); X2 = X1+y_kspace-1;
cropR = floor(range(positions(:,2))/pixel_size);
cropC = floor(range(positions(:,1))/pixel_size);
%% create initial aperture?and object guesses
if aperture == 0
    aperture = single(((2*makeCircleMask(round(aperture_radius./pixel_size),little_area))));
    if strcmp(ApType,'F')
        aperture = my_ifft(aperture);
    end
    initial_aperture = aperture;
else
    aperture = single(aperture);
    initial_aperture = aperture;
end
if big_obj == 0
    %big_obj = single(ones(bigx,bigy)).*exp(1i*(rand(bigx,bigy)));
    big_obj = obj_scale*ones(bigx,bigy,'single');
    %big_obj = zeros(bigx,bigy);
    initial_obj = big_obj;
else
    big_obj = single(big_obj);
    initial_obj = big_obj;
end

%{
[XX,YY] = meshgrid(1:bigx,1:bigy);
X_cen = floor(bigx/2); Y_cen = floor(bigy/2);
R2 = (XX-X_cen).^2 + (YY-Y_cen).^2;
Kfilter = exp(-R2/(2*400)^2);
Kfilter = Kfilter/max(Kfilter(:));
%}

fourier_error = zeros(iterations,nApert);
sum_dp = zeros(nApert,1);
if gpu
    display('========DR reconstructing with GPU========')
    diffpats = gpuArray(diffpats);
    fourier_error = gpuArray(fourier_error);
    big_obj = gpuArray(big_obj);
    aperture = gpuArray(aperture);
    sum_dp = gpuArray(sum_dp);
else
    display('========DR reconstructing with CPU========')
end


for aper=1:nApert
    diffpat_i = diffpats(:,:,aper);
    missing_data = diffpat_i == -1;
    sum_dp(aper) = sum(sum(diffpat_i(~missing_data)));
end

best_err = 100; % check to make sure saving reconstruction with best error
Z = diffpats.*exp(rand(N1,N2,nApert));

%weights = init_weight + 0.5*(final_weight-init_weight)*round(2*((1:iterations)/iterations).^order);
weights = init_weight + (final_weight-init_weight)*((1:iterations)/iterations).^order;
alpha_temp = alpha;
alpha=0;

%% Main ePIE itteration loop
disp('========beginning reconstruction=======');
tic
for t = 1:iterations
    %weight = weights(t);
    weight = 0.1;
    %if t==20, alpha = alpha_temp; end
    if t==20, alpha = 0.4+ (alpha_temp-0.4)*(iterations-t)/iterations; end
    
    %% sequential
    for aper = randperm(nApert)
        %bigObjShifted = circshift(big_obj, [-1*(centrey(aper) - centBig) -1*(centrex(aper) - centBig)]);
        %u_old = croppedOut(bigObjShifted,y_kspace); % size 384x384
        u_old = big_obj(Y1(aper):Y2(aper), X1(aper):X2(aper));
        %object_max = max(abs(u_old(:))).^2;
        probe_max = max(abs(aperture(:))).^2;
        %% Create new exitwave
        Pu_old = u_old.*(aperture);
        z_u = fft2(Pu_old);
        check_dp = abs(z_u);
        
        z = Z(:,:,aper);
        z_F = (1+alpha)*z_u - alpha*z;
        
        diffpat_i = diffpats(:,:,aper);
        missing_data = diffpat_i == -1;
        k_fill = z_F(missing_data);
        
        % update z: tunning weight
        phase = angle(z_F);
        z_F =  (1-weight)*exp(1i*phase) .* diffpats(:,:,aper) + weight*z_F ;
        z_F(missing_data) = k_fill;
        z = z_F + alpha*(z- z_u);
        Z(:,:,aper)=z;
        
        %% Update the object
        Pu_new = ifft2(z);
        diff = Pu_old - Pu_new;
        dt = beta_obj/probe_max;
        %u_temp = (u_old + dt*Pu_new.*conj(aperture)) ./ ( 1 + dt*abs(aperture).^2 );
        u_temp = ( ((1-beta_obj)/dt)*u_old + Pu_new.*conj(aperture)) ./ ( (1-beta_obj)/dt + abs(aperture).^2 );
        
        if do_posi, u_new = max(0,real(u_new));
        else u_new = u_temp; end
        
        big_obj(Y1(aper):Y2(aper), X1(aper):X2(aper)) = u_new;
        %% Update the probe
        if t > 20
        if update_aperture && t > iterations-freeze_aperture
            object_max = max(abs(u_new(:))).^2;
            new_beta_ap = beta_ap*sqrt((iterations-t+1)/iterations);
            ds = new_beta_ap./(object_max+eps);
            
            if semi_implicit_P
                aperture = ((1-new_beta_ap)*aperture + ds*Pu_new.*conj(u_old)) ./ ( (1-new_beta_ap) + ds*abs(u_old).^2 );
            else aperture = aperture - ds*conj(u_old).*(diff);
            end
            
            switch probe_mask_flag
                case 0
                case 1 %real space
                    aperture = aperture .* R_probe_mask;
                case 2 %k-space
                    aperture = fftshift(ifft2(fft2(ifftshift(aperture)) .* F_probe_mask));
                otherwise
                    error('probe mask flag can only be 0,1 or 2');
            end
        end
        
    end
        
        if sum_dp(aper)
            fourier_error(t,aper) = sum(abs( diffpat_i(~missing_data) - check_dp(~missing_data) )) ./ sum_dp(aper);
        end
    end    
    
    if probeNorm == 1
        scale = max(max(abs(aperture)));
        aperture = aperture/scale;
    end
    %{
    if t>iterations-100
    z = fftshift(fft2(big_obj)).*Kfilter;
    big_obj = ifft2(ifftshift(z));
    end
    %}
    
    %% show results
    if  mod(t,showim) == 0 && imout == 1;
        %         figure(333)
        %         %hsv_big_obj = make_hsv(big_obj,1);
        %         hsv_aper = make_hsv(aperture,1);
        %         subplot(2,2,1)
        %         imagesc(abs(big_obj)); axis image; colormap gray; title(['reconstruction pixel size = ' num2str(pixel_size)] )
        %         subplot(2,2,2)
        %         imagesc(hsv_aper); axis image; colormap gray; title('aperture single'); colorbar
        %         subplot(2,2,3)
        %         plot(errors); ylim([0,0.5]);
        %         subplot(2,2,4)
        %         imagesc(log(fftshift(check_dp))); axis image
        %         drawnow
        %figure(39);
        %img_i = abs(big_obj(276:555,276:555));
        %images(:,:,t) = img_i;
        %imagesc(img_i); axis image; colormap(jet); colorbar;
        %drawnow;
        %drawnow
        %writeVideo(writerObj,getframe(h)); % take screen shot
        
        hsv_ap = make_hsv(aperture,1);
        %figure(33); img(croppedOut(big_obj,[cropR,cropC]), ['DR ',num2str(t)], 'colormap', 'gray');
        %figure(33); img(big_obj, ['DR ',num2str(t)], 'colormap', 'gray');
        figure(33); imagesc(abs(big_obj)), title(['DR ',num2str(t)]); colormap gray; axis image
        figure(34); imagesc(hsv_ap); title(num2str(t)); colormap gray; axis image;
        drawnow
    end
    
    mean_err = sum(fourier_error(t,:),2)/nApert;
    if best_err > mean_err
        best_obj = big_obj;
        best_err = mean_err;
    end
    if save_intermediate == 1 && mod(t,save_inter_freq) == 0
        save(sprintf('%s%s_itt%04d.mat',save_string,filename,t),...%[save_string 'inter_output_rPIE_ID_',jobID,'_itt_',num2str(itt),'.mat'],...
            'best_obj','big_obj','aperture','fourier_error','-v7.3');
    end
    fprintf('%d. Error = %f\n',t,mean_err);
    if t>20
            [~,~,err_amp,err_angle] = quancomp(big_obj,object);
    hist.amp = [hist.amp err_amp];
    hist.angle = [hist.angle err_angle];
    hist.Rfactor = [hist.Rfactor mean_err];
    end
end
toc
disp('======reconstruction finished=======')

if gpu
    fourier_error = gather(fourier_error);
    best_obj = gather(best_obj);
    aperture = gather(aperture);
end

if saveOutput
    save(sprintf('%soutput_%s.mat',save_string,filename),...%[save_string 'best_obj_' filename '.mat'],'best_obj',
        'big_obj','best_obj','aperture','fourier_error');
end

%% Function for converting positions from experimental geometry to pixel geometry

%     function [positions] = convert_to_pixel_positions(positions,pixel_size)
%         positions = positions./pixel_size;
%         positions(:,1) = (positions(:,1)-min(positions(:,1)));
%         positions(:,2) = (positions(:,2)-min(positions(:,2)));
%         positions(:,1) = (positions(:,1)-round(max(positions(:,1))/2));
%         positions(:,2) = (positions(:,2)-round(max(positions(:,2))/2));
%         positions = round(positions);
%         bigx =little_area + max(positions(:))*2+10; % Field of view for full object
%         bigy = little_area + max(positions(:))*2+10;
%         big_cent = floor(bigx/2)+1;
%         positions = positions+big_cent;
%
%
%     end
    function [shiftDistancesX, shiftDistancesY, truePositions, positions, bigx, bigy] = ...
            convert_to_pixel_positions_testing3(positions,pixel_size,little_area)
        
        
        positions = positions./pixel_size;
        positions(:,1) = (positions(:,1)-min(positions(:,1)));
        positions(:,2) = (positions(:,2)-min(positions(:,2)));
        truePositions(:,1) = (positions(:,1) - max(positions(:,1))/2);
        truePositions(:,2) = (positions(:,2) - max(positions(:,2))/2);
        positions(:,1) = (positions(:,1)-round(max(positions(:,1))/2));
        positions(:,2) = (positions(:,2)-round(max(positions(:,2))/2));
        
        positions = round(positions);
        bigx =little_area + max(positions(:))*2+10; % Field of view for full object
        bigy = little_area + max(positions(:))*2+10;
        %         bigx = little_area;
        %         bigy = little_area;
        big_cent = floor(bigx/2)+1;
        positions = positions+big_cent;
        truePositions = truePositions + big_cent;
        
        shiftDistancesX = truePositions(:,1) - positions(:,1);
        shiftDistancesY = truePositions(:,2) - positions(:,2);
        
        
        
    end

    function [shiftDistancesX, shiftDistancesY, truePositions, positions, bigx, bigy] = ...
            convert_to_pixel_positions_testing4(positions,pixel_size,little_area)
        
        
        positions = positions./pixel_size;
        positions(:,1) = (positions(:,1)-min(positions(:,1)));
        positions(:,2) = (positions(:,2)-min(positions(:,2)));
        truePositions(:,1) = (positions(:,1) - round(max(positions(:,1))/2));
        truePositions(:,2) = (positions(:,2) - round(max(positions(:,2))/2));
        
        positions = round(truePositions);
        bigx =little_area + max(positions(:))*2+10; % Field of view for full object
        bigy = little_area + max(positions(:))*2+10;
        %         bigx = little_area;
        %         bigy = little_area;
        big_cent = floor(bigx/2)+1;
        positions = positions+big_cent;
        truePositions = truePositions + big_cent;
        
        shiftDistancesX = truePositions(:,1) - positions(:,1);
        shiftDistancesY = truePositions(:,2) - positions(:,2);
        
        
        
    end

    function [pixelPositions, bigx, bigy] = ...
            convert_to_pixel_positions_testing5(positions,pixel_size,little_area)
        
        
        pixelPositions = positions./pixel_size;
        pixelPositions(:,1) = (pixelPositions(:,1)-min(pixelPositions(:,1))); %x goes from 0 to max
        pixelPositions(:,2) = (pixelPositions(:,2)-min(pixelPositions(:,2))); %y goes from 0 to max
        pixelPositions(:,1) = (pixelPositions(:,1) - round(max(pixelPositions(:,1))/2)); %x is centrosymmetric around 0
        pixelPositions(:,2) = (pixelPositions(:,2) - round(max(pixelPositions(:,2))/2)); %y is centrosymmetric around 0
        
        bigx = little_area + round(max(pixelPositions(:)))*2+10; % Field of view for full object
        bigy = little_area + round(max(pixelPositions(:)))*2+10;
        
        big_cent = floor(bigx/2)+1;
        
        pixelPositions = pixelPositions + big_cent;
        
        
    end

%% 2D gaussian smoothing of an image

    function [smoothImg,cutoffRad]= smooth2d(img,resolutionCutoff)
        
        Rsize = size(img,1);
        Csize = size(img,2);
        Rcenter = round((Rsize+1)/2);
        Ccenter = round((Csize+1)/2);
        a=1:1:Rsize;
        b=1:1:Csize;
        [bb,aa]=meshgrid(b,a);
        sigma=(Rsize*resolutionCutoff)/(2*sqrt(2));
        kfilter=exp( -( ( ((sqrt((aa-Rcenter).^2+(bb-Ccenter).^2)).^2) ) ./ (2* sigma.^2) ));
        kfilter=kfilter/max(max(kfilter));
        kbinned = my_fft(img);
        
        kbinned = kbinned.*kfilter;
        smoothImg = my_ifft(kbinned);
        
        [Y, X] = ind2sub(size(img),find(kfilter<(exp(-1))));
        
        Y = Y-(size(img,2)/2);
        X = X-(size(img,2)/2);
        R = sqrt(Y.^2+X.^2);
        cutoffRad = ceil(min(abs(R)));
    end

%% Fresnel propogation
    function U = fresnel_advance (U0, dx, dy, z, lambda)
        % The function receives a field U0 at wavelength lambda
        % and returns the field U after distance z, using the Fresnel
        % approximation. dx, dy, are spatial resolution.
        
        k=2*pi/lambda;
        [ny, nx] = size(U0);
        
        Lx = dx * nx;
        Ly = dy * ny;
        
        dfx = 1./Lx;
        dfy = 1./Ly;
        
        u = ones(nx,1)*((1:nx)-nx/2)*dfx;
        v = ((1:ny)-ny/2)'*ones(1,ny)*dfy;
        
        O = my_fft(U0);
        
        H = exp(1i*k*z).*exp(-1i*pi*lambda*z*(u.^2+v.^2));
        
        U = my_ifft(O.*H);
    end

%% Make a circle of defined radius

    function out = makeCircleMask(radius,imgSize)
        
        
        nc = imgSize/2+1;
        n2 = nc-1;
        [xx, yy] = meshgrid(-n2:n2-1,-n2:n2-1);
        R = sqrt(xx.^2 + yy.^2);
        out = R<=radius;
    end

%% Function for creating HSV display objects for showing phase and magnitude
%  of a reconstruction simaultaneously

    function [hsv_obj] = make_hsv(initial_obj, factor)
        
        [sizey,sizex] = size(initial_obj);
        hue = angle(initial_obj);
        
        value = abs(initial_obj);
        hue = hue - min(hue(:));
        hue = (hue./max(hue(:)));
        value = (value./max(value(:))).*factor;
        hsv_obj(:,:,1) = hue;
        hsv_obj(:,:,3) = value;
        hsv_obj(:,:,2) = ones(sizey,sizex);
        hsv_obj = hsv2rgb(hsv_obj);
    end
%% Function for defining a specific region of an image

    function [roi] = get_roi(image, centrex,centrey,crop_size)
        
        bigy = size(image,1);
        bigx = size(image,2);
        
        half_crop_size = floor(crop_size/2);
        if mod(crop_size,2) == 0
            roi = {centrex - half_crop_size:centrex + (half_crop_size - 1);...
                centrey - half_crop_size:centrey + (half_crop_size - 1)};
            
        else
            roi = {centrex - half_crop_size:centrex + (half_crop_size);...
                centrey - half_crop_size:centrey + (half_crop_size)};
            
        end
    end

%% Fast Fourier transform function
    function kout = my_fft(img)
        kout = fftshift(fftn((ifftshift(img))));
    end
%% Inverse Fast Fourier transform function
    function realout = my_ifft(k)
        realout =fftshift((ifftn(ifftshift(k))));
        %realout =ifftshift((ifftn(fftshift(k))));
        
    end

    function ROI = croppedOut(largeArray,cropSize)
        %     n = size(largeArray,1);
        %     nc = round((n+1)/2);
        %
        %     cropVec = 1:cropSize;
        %     cropC = round((cropSize+1)/2);
        %     cropVec = cropVec - cropC + nc;
        %     ROI = largeArray(cropVec, cropVec);
        n = size(largeArray);
        nc = round((n+1)/2);
        
        if length(cropSize) == 1
            cropSize = repmat(cropSize,length(n),1);
        end
        for i = 1:length(n)
            vec = 1:cropSize(i);
            cropc = round((cropSize(i)+1)/2);
            cropVec{i} = single(vec - cropc + nc(i));
        end
        
        if length(n) == 2
            ROI = largeArray(cropVec{1}, cropVec{2});
        elseif length(n) == 3
            ROI = largeArray(cropVec{1}, cropVec{2}, cropVec{3});
        end
    end
end


