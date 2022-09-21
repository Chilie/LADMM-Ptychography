addpath('functions\');
addpath('data\');
ap_radius    = 50;%50; % 50
scaning_steps = [21,25,27]; % 51
for i_s = 1:3
    scaning_step = scaning_steps(i_s);
scan_type = 'spiral';%'spiral';%'grid'; % 'spiral';
sigma = 100; %30;
%phases = importdata('model.mat'); small_phase = phases(129:256,129:256);
phase = im2double(imread('pepper.png')); phase = double(phase(:,:,1));  phase = padarray(phase,[128,128],0,'both');
model = im2double(imread('cameraman.png')); model = double(model(:,:,1)); model = padarray(model,[128,128],0,'both');

% [a, ~, centerx, centery] = make_apertures(model,scaning_step,ap_radius,1e6,'grid',3,120,120); % 3
% [a, ~, centerx, centery] = make_apertures(model,scaning_step,ap_radius,sigma,'grid',0,120,120); % 3
[a, ~, centerx, centery] = make_apertures(model,scaning_step,ap_radius,100,'spiral',0,120,120);
[N1,N2,nProbes] = size(a);
%a = a.*repmat(exp(1i*(rand(N1,N2)-.5)),[1,1,nProbes]);

%% generate diffraction patterns with poisson noise
% flux
flux=1e8;
dp = zeros([N1,N2,nProbes]);
dp0 = zeros([N1,N2,nProbes]);
for ii = 1:nProbes
    dp0(:,:,ii) = abs(fftshift(fftn(model.*exp(1i*(2*pi*phase-pi)).*a(:,:,ii)))).^2;
%     dpi = dp0(:,:,ii);
%     num_sum = sum(dpi(:));
%     if num_sum
%     scale = flux/sum(dpi(:));
%     dp(:,:,ii) = poissrnd(dpi*scale)./scale;
%     else
%     dp(:,:,ii) = dp0(:,:,ii);
%     end
    dp(:,:,ii) = dp0(:,:,ii);
end
%% inputs
object = model.*exp(1i*(2*pi*phase-pi)); object = object(129:256,129:256);

save_path = 'Pattern_images'; % '_noiseless'
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

%% plot the scanning pattern
positions = [centerx' centery'];
[m,n] = size(dp(:,:,1));
Scan_pattern = plot_the_scanning_pattern(positions, m,n);
Inner = imgaussfilt(Scan_pattern(128:end-128,128:end-128),3);
figure;imagesc(Inner);
Inner = uint8(255 * mat2gray(Inner));
imwrite(Inner, fullfile(save_path, [scan_type '_scan_step_' num2str(scaning_step) '_pattern.png']));
temp = a(:,:,1);

probe_gt = temp(centerx(1)-123:centerx(1)+132, centery(1)-123:centery(1)+132);
Inner = uint8(255 * mat2gray(probe_gt));
imwrite(ind2rgb(Inner,parula(256)), fullfile(save_path, [scan_type '_scan_step_' num2str(scaning_step) '_pattern_probe.png']));
end
%% normlized object