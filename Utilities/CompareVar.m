Var = Var_SigP2;

Object = 'Cell';

NumObj = size(Var.Measurements.(Object).PixelList,2);
ObjImg = zeros(Var.Analysis.ImgSize);
Tpt = 10;


for obj = 1:NumObj
    ObjImg(Var.Measurements.(Object).PixelList{Tpt, obj}) = obj;
end

figure(1)
imagesc(ObjImg)

ObjImgNorm = zeros(Var.Analysis.ImgSize);
ObjImgNorm(ObjImg>0) = .8;

%% ___________________________

Var = Var_AP2;

Object = 'Cell';

NumObj = size(Var.Measurements.(Object).PixelList,2);
ObjImg_2 = zeros(Var.Analysis.ImgSize);
%Tpt = 1;


for obj = 1:NumObj
    ObjImg_2(Var.Measurements.(Object).PixelList{Tpt, obj}) = obj;
end

figure(2)
imagesc(ObjImg_2)



ObjImgNorm_2 = zeros(Var.Analysis.ImgSize);
ObjImgNorm_2(ObjImg_2>0) = 0.8;


RGBimg(:,:,3) = ObjImgNorm_2;
RGBimg(:,:,1) = ObjImgNorm;
figure(3)
image(RGBimg)

