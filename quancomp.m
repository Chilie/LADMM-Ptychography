function [norm_amp,norm_angle,err_amp, err_angle] = quancomp(big_obj,object)
%UNTITLED10 compute the metrics
%   此处显示详细说明
[m,n] = size(object);
correlation1 = normxcorr2(abs(object),abs(big_obj));
h1 = round(size(big_obj)/2);
max1 = max(max(abs(correlation1(h1-m:h1+m-1,h1-n:h1+n-1)) ));
I = find(abs(correlation1)==max1);
[I1,I2] = ind2sub(size(correlation1),I);

object1 = big_obj(I1-size(object,1)+1:I1, I2-size(object,2)+1:I2 );
%object1 = big_obj(I1-size(object,1)+2:I1+1, I2-size(object,2)+2:I2+1 );

shift1 = norm(object1,'fro');
norm_amp = abs(object1/shift1);
err_amp = norm(norm_amp - abs(object),'fro');
% shift1 = sum(conj(object1(:)).*object(:)); shift1 = shift1/norm(shift1);
norm_angle = angle(object1/shift1); 
err_angle = norm(norm_angle - angle(object),'fro')/norm(angle(object),'fro');
end

