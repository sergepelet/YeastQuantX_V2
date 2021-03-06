
mkdir YQ_VX/src
mkdir YQ_VX/distrib

mcc -o YQ_VX -W main:YQ_VX -T link:exe -d /users/spelet1/YeastQuant/YQ_VX/src...
    -N -p images -p signal -p stats -w enable:specified_file_mismatch -w enable:repeated_file -w enable:switch_ignored -w enable:missing_lib_sentinel -w enable:demo_license -R -nodisplay -R -singleCompThread ...
    -v /users/spelet1/YeastQuant/WrapFrameAnalysis.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/CalculatedObject.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/CombineData.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/CombineFrame.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/ExpandObject.m ...
    -a /users/spelet1/YeastQuant/FrameAnalysis.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/Img2MovieRGB.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/LoadFrameMeas.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/LoadImage.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/MeasureObject.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/ObjectTransfer.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/RemoveMeasurements.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/ResolveNuclAssignment.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/SecondaryObject.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/SegmentCellX.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/SegmentNuclBF.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/SegmentYeastFluo.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/SegmentYeastFluoMedia.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/SegmentYeastPhase.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/SplitGroupHough.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/SubtractBackground.m ...
    -a /users/spelet1/YeastQuant/CommonProcess/TrackObjectsSingle.m ...
      -a /users/spelet1/YeastQuant/CellX/core/mclasses/CellXConfiguration.m ...
        -a /users/spelet1/YeastQuant/CellX/core/mclasses/CellXFileSet.m ...
  -a /users/spelet1/YeastQuant/CellX/core/mclasses/CellXMembraneDetector.m ...
  -a /users/spelet1/YeastQuant/CellX/core/mclasses/CellXSeed.m ...
  -a /users/spelet1/YeastQuant/CellX/core/mclasses/CellXSeedIntersectionResolver.m ...
  -a /users/spelet1/YeastQuant/CellX/core/mclasses/CellXSeedValidator.m ...
    -a /users/spelet1/YeastQuant/CellX/core/mclasses/CellXSegmenter.m ...
  -a /users/spelet1/YeastQuant/CellX/core/mex/computeRayConvolution.mexa64 ...
  -a /users/spelet1/YeastQuant/CellX/core/mex/edgeJoin.mexa64 ...
  -a /users/spelet1/YeastQuant/CellX/core/mex/findRegionCentersInCircularArray.mexa64 ...
  -a /users/spelet1/YeastQuant/CellX/core/mex/initGridGraph.mexa64 ...
  -a /users/spelet1/YeastQuant/CellX/core/mex/maxflow/maxflow.m ...
  -a /users/spelet1/YeastQuant/CellX/core/mex/maxflow/maxflowmex.mexa64 ...
  -a /users/spelet1/YeastQuant/CellX/core/mfunctions/CircularHough_Grd.m
  
    