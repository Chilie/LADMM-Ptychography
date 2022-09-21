addpath('functions');
addpath('data');
ap_radius    = 50;%50; % 50
scaning_step = 41; % 51
scan_type = 'grid';%'spiral';%'grid'; % 'spiral';
sigma = 100; %30;
%phases = importdata('model.mat'); small_phase = phases(129:256,129:256);
phase = im2double(imread('pepper.png')); phase = double(phase(:,:,1));  phase = padarray(phase,[128,128],0,'both');
model = im2double(imread('cameraman.png')); model = double(model(:,:,1)); model = padarray(model,[128,128],0,'both');

% [a, ~, centerx, centery] = make_apertures(model,scaning_step,ap_radius,1e6,'grid',3,120,120); % 3
[a, ~, centerx, centery] = make_apertures(model,scaning_step,ap_radius,sigma,'grid',0,120,120); % 3
% [a, ~, centerx, centery] = make_apertures(model,scaning_step,ap_radius,100,'spiral',0,120,120);
[N1,N2,nProbes] = size(a);
%a = a.*repmat(exp(1i*(rand(N1,N2)-.5)),[1,1,nProbes]);

 set(gcf,'Visible','off');              
 set(0,'DefaultFigureVisible','off');

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


%% plot the scanning pattern
positions = [centerx' centery'];
[m,n] = size(dp(:,:,1));
Scan_pattern = plot_the_scanning_pattern(positions, m,n);
Inner = imgaussfilt(Scan_pattern(128:end-128,128:end-128),3);
figure;imagesc(Inner);
Inner = uint8(255 * mat2gray(Inner));
% imwrite(Inner, fullfile(save_path, 'pattern.png'));
temp = a(:,:,1);
% probe_gt = temp(centerx(1)-123:centerx(1)+132, centery(1)-123:centery(1)+132);
% Inner = uint8(255 * mat2gray(probe_gt));
% imwrite(ind2rgb(Inner,parula(256)), fullfile(save_path, 'pattern_probe.png'));
%% normlized object

Test_No = 20;
for I_test = 1:Test_No
save_path = ['Fixed_init_Revised_tmp_probe_nongau_sigma_' num2str(sigma) '_' scan_type '_ap_radius_' num2str(ap_radius) '_scan_step_' num2str(scaning_step) '_noiseless/Run_' num2str(I_test)]; % '_noiseless'
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

object_norm = object/norm(object,'fro');
ePIE_inputs.GpuFlag = 0;
ePIE_inputs.Patterns = dp;
ePIE_inputs.Positions = [centerx' centery'];
ePIE_inputs.FileName = 'ePIE_cameraman';
ePIE_inputs.PixelSize = 1;
ePIE_inputs.InitialObj = 0;
ePIE_inputs.ApRadius = ap_radius; %radius of aperture
ePIE_inputs.ApType  = 'R';        %aperture type: R: real unique circle, F: IFFT of a unique circle
ePIE_inputs.InitialAp = 0;
ePIE_inputs.Iterations = 220;
ePIE_inputs.showim = 0;
ePIE_inputs.updateAp = 1;
ePIE_inputs.object = object_norm;

noise = sum(sum(abs(sqrt(dp(:)) - sqrt(dp0(:))))) / sum(sum(sqrt(dp0(:))));
fprintf('Noise = %f\n',noise);

%% obtain the size of bigx bigy
[bigx, bigy] = determinesize(ePIE_inputs);
ePIE_inputs.InitialObj = single(rand(bigx,bigy)).*exp(1i*(rand(bigx,bigy)));

%% ePIE reconstruction
ePIE_inputs.Iterations = 220;
ePIE_inputs.do_posi = 0;
ePIE_inputs.FileName = 'ePIE_test';
ePIE_inputs.Verbosity = 1;
[big_obj_0, aperture_0, fourier_error, initial_obj, initial_aperture, hist_epie]  = ePIE(ePIE_inputs,1,0.01);
% [big_obj, aperture, fourier_error, initial_obj, initial_aperture,hist]  = LADMM_PIE(ePIE_inputs,1,0.01);
% [big_obj, aperture, fourier_error, initial_obj, initial_aperture]  = ADMM_PIE(ePIE_inputs,1,0.01);
ePIE_inputs.FileName = 'LADMM_test';
[big_obj_ladmm, aperture_ladmm, fourier_error, initial_obj, initial_aperture,hist_ladmm]  = Adaptive_LADMM_PIE(ePIE_inputs,1,0.01);

ePIE_inputs.FileName = 'ADMM_test';
[big_obj_admm, aperture_admm, fourier_error, initial_obj, initial_aperture,hist_admm]  = ADMM_PIE_mutiple(ePIE_inputs,1,0.01);
% 
ePIE_inputs.FileName = 'PhiBE_test';
[big_obj_phibe, aperture_phibe, fourier_error, initial_obj, initial_aperture,hist_phibe]  = PhiBE(ePIE_inputs,1,0.01);
% [big_obj, aperture, fourier_error, initial_obj, initial_aperture]  = ADMM_PIE(ePIE_inputs,1,0.01);
% rPIE reconstruction
ePIE_inputs.Iterations = 220;
ePIE_inputs.FileName = 'rPIE_test';
[big_obj2,aperture2,fourier_error2,initial_obj2,initial_aperture2,hist_rpie] = rPIE(ePIE_inputs,0.1,0.02);

% DR reconstruction
ePIE_inputs.Iterations = 220;
ePIE_inputs.FileName = 'DR_test';
[big_obj3,aperture3,fourier_error3,initial_obj3,initial_aperture3,hist_dr] = DRb(ePIE_inputs,0.7,0.01,0.9);

Tabs = {'big_obj_0','big_obj_ladmm','big_obj_admm','big_obj_phibe'};
Tabs_probe = {'aperture_0','aperture_ladmm','aperture_admm','aperture_phibe'};
for ik = 1:4
eval(['big_obj=' Tabs{ik} ';']);
eval(['aperture=' Tabs_probe{ik} ';']);
%% result of ePIE
correlation1 = normxcorr2(abs(object),abs(big_obj));
h1 = round(size(big_obj)/2);
% max1 = max(max(abs(correlation1(h1-128:h1+127,h1-128:h1+127)) ));
max1 = max(max(abs(correlation1(h1-200:h1+199,h1-200:h1+199)) ));
I = find(abs(correlation1)==max1);
[I1,I2] = ind2sub(size(correlation1),I);

object1 = big_obj(I1-size(object,1)+1:I1, I2-size(object,2)+1:I2);
%object1 = big_obj(I1-size(object,1)+2:I1+1, I2-size(object,2)+2:I2+1 );

figure(11); imagesc(abs(object1)); axis image; colormap gray;
set(gca, 'visible', 'off');
Inta = uint8(255 * mat2gray(abs(object1)));

shift1 = sum(conj(object1(:)).*object(:)); shift1 = shift1/norm(shift1);
angle1 = angle(object1*shift1); 
figure(12); imagesc(angle1); axis image; colormap gray;
set(gca, 'visible', 'off');
Phas = uint8(255 * mat2gray(angle1));
imwrite(Inta, fullfile(save_path, [Tabs{ik} '_amp.png']));
imwrite(Phas, fullfile(save_path, [Tabs{ik} '_ang.png']));

rec_probe = uint8(255 * mat2gray(abs(aperture)));
imwrite(ind2rgb(rec_probe,parula(256)), fullfile(save_path, [Tabs_probe{ik} '_prob.png']));
end

%% result of rPIE
correlation2 = normxcorr2(abs(object),abs(big_obj2));
max1 = max(max(abs(correlation2(h1-128:h1+127,h1-128:h1+127)) ));
I = find(abs(abs(correlation2)-max1)<eps);
[I1,I2] = ind2sub(size(correlation2),I);

object2 = big_obj2(I1-size(object,1)+1:I1, I2-size(object,2)+1:I2 );
%object2 = big_obj2(I1-size(object,1)+2:I1+1, I2-size(object,2)+2:I2+1 );


figure(21); imagesc(abs(object2)); axis image; colormap gray;
set(gca, 'visible', 'off');
Inta = uint8(255 * mat2gray(abs(object2)));
shift2 = sum(conj(object2(:)).*object(:)); shift2 = shift2/norm(shift2);
angle2 = angle(object2*shift2); 
figure(22); imagesc(angle2); axis image; colormap gray;
set(gca, 'visible', 'off');
Phas = uint8(255 * mat2gray(angle2));
imwrite(Inta, fullfile(save_path, ['rpie' '_amp.png']));
imwrite(Phas, fullfile(save_path, ['rpie' '_ang.png']));

rec_probe = uint8(255 * mat2gray(abs(aperture2)));
imwrite(ind2rgb(rec_probe,parula(256)), fullfile(save_path, 'rpie_prob.png'));
%% result of sDR
correlation3 = normxcorr2(abs(object),abs(big_obj3));
max1 = max(max(abs(correlation3(h1-128:h1+127,h1-128:h1+127)) ));
I = find(abs(correlation3)==max1);
[I1,I2] = ind2sub(size(correlation3),I);

object3 = big_obj3(I1-size(object,1)+1:I1, I2-size(object,2)+1:I2 );
%object3 = big_obj3(I1-size(object,1)+2:I1+1, I2-size(object,2)+2:I2+1 );


figure(31); imagesc(abs(object3)); axis image; colormap gray;
set(gca, 'visible', 'off');
Inta = uint8(255 * mat2gray(abs(object3)));
shift3 = sum(conj(object3(:)).*object(:)); shift3 = shift3/norm(shift3);
angle3 = angle(object3*shift3); 
figure(32); imagesc(angle3); axis image; colormap gray;
set(gca, 'visible', 'off');
Phas = uint8(255 * mat2gray(angle3));

imwrite(Inta, fullfile(save_path, ['sdr' '_amp.png']));
imwrite(Phas, fullfile(save_path, ['sdr' '_ang.png']));
rec_probe = uint8(255 * mat2gray(abs(aperture3)));
imwrite(ind2rgb(rec_probe,parula(256)), fullfile(save_path, 'sdr_prob.png'));
%%
sphase = im2double(imread('pepper.png')); sphase = double(sphase(:,:,1)); 
smodel = im2double(imread('cameraman.png')); smodel = double(smodel(:,:,1));
figure(101); imagesc(sphase); colormap gray;axis image;
set(gca, 'visible', 'off');
figure(102); imagesc(smodel); colormap gray;axis image;
set(gca, 'visible', 'off');

imwrite(uint8(255 * mat2gray(smodel)), fullfile(save_path, ['gt' '_amp.png']));
imwrite(uint8(255 * mat2gray(sphase)), fullfile(save_path, ['gt' '_ang.png']));

save(fullfile(save_path, 'data.mat'), 'hist_*','big_obj*','aperture*');
end
