function Big_m = plot_the_scanning_pattern(positions, m,n)
% input [m,n] the size of the difffraction pattern
y_kspace = [m,n];
pixel_size = 1;
little_area = y_kspace; % Region of reconstruction for placing back into full size image

%% Get centre positions for cropping (should be a 2 by n vector)
%positions = convert_to_pixel_positions(positions,pixel_size);
[pixelPositions, bigx, bigy] = ...
    convert_to_pixel_positions_testing5(positions,pixel_size,little_area);
centrey = round(pixelPositions(:,2));
centrex = round(pixelPositions(:,1));
centBig = round((bigx+1)/2);
Y1 = centrey - floor(y_kspace/2); Y2 = Y1+y_kspace-1;
X1 = centrex - floor(y_kspace/2); X2 = X1+y_kspace-1;

Num = length(Y1);
Big_m = zeros(bigx);

for kk = 1:Num
   Big_m(centrey(kk),centrex(kk)) = 1; 
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

end