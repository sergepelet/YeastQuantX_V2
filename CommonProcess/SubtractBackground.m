function Var = SubtractBackground(Var, CallNum)
tic
if nargin == 1
    CallNum = 1;
end
%Assigne Ref and input images
Image = Var.Img.(Var.Analysis.SubtractIn{CallNum});
BackgroundImage = Var.Img.(Var.Analysis.SubtractRef{CallNum});

%Grow background image
if Var.Analysis.SubtractGrow{CallNum} ~= 0
    SE = strel('disk',Var.Analysis.SubtractGrow{CallNum});
    BackgroundImage = imdilate(BackgroundImage,SE);
end

    %Generate histogram of pixel intensities outside of the objects
    %Remove 2%  hi and low pixels
    Xhist = linspace(prctile(Image(BackgroundImage == 0),0.2),prctile(Image(BackgroundImage == 0),99.8),500);
    [Yhist] = hist(Image(BackgroundImage == 0),Xhist);
try
    %Fit with gaussian
    % generate guess for start of fit
%     [MaxInt, MaxInd] = max(Yhist);
%     MedImage = median(Image(BackgroundImage == 0));
%     coef0 = [MaxInt, MedImage , 10];
    FitOpts = fitoptions('method','NonlinearLeastSquares'); %,'StartPoint',coef0);
    FitType = fittype('gauss1');
    FitResult = fit(Xhist', Yhist', FitType,FitOpts);
    YFit = FitResult(Xhist);
    CoefFit = coeffvalues(FitResult);
    MeanGauss = CoefFit(2);
catch
    [MaxInt, MaxInd] = max(Yhist);
    MeanGauss = Xhist(MaxInd);
end
    

%Remove Background intensity from image
RemainingImage = Image - MeanGauss;


%Set Image to zero in the background
%BackgroundIntensity = mean(Image(BackgroundImage == 0))
%RemainingImage = Image - BackgroundIntensity;

%%% Display %%%
if strcmp(Var.Figure.Display, 'on')
    FigNum = find(strcmp(Var.Figure.List, 'SubBackground'));
    figure(FigNum(CallNum))
    subplot(2,2,1); imagesc(Image);title('Intensity Image')
    subplot(2,2,2); imagesc(BackgroundImage); title('Reference Image');
    subplot(2,2,4); imagesc(RemainingImage); title('Background Subtracted Image');
    subplot(2,2,3);
    plot(Xhist, Yhist)
    if exist('YFit', 'var')
        hold on
        plot(Xhist, YFit, 'r')
        hold off
    end
end


%%% save %%%
Var.Img.(Var.Analysis.SubtractOut{CallNum}) = RemainingImage;

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;