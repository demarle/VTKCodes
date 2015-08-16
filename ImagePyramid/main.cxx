/*=========================================================================

  Program:   Visualization Toolkit
  Module:    main.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//Program to demonstrate Eulerian Motion Magnification
//by Ashray Malhotra
//for more information see http://www.kitware.com/blog/home/post/952

#define YChannel 0
#define IChannel 1
#define QChannel 2

// #define debug

#include "vtkImagePyramid.h"

#include <vtkGlobFileNames.h>
#include <vtkImageAppendComponents.h>
#include <vtkImageConvolve.h>

#include <vtkImageData.h>
#include <vtkImageDifference.h>

#include <vtkImageExtractComponents.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageMathematics.h>
#include <vtkImageReader2.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageResize.h>
#include <vtkImageRGBToYIQ.h>
#include <vtkImageWeightedSum.h>
#include <vtkImageYIQToRGB.h>
#include <vtkPNGWriter.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

#include <vtksys/SystemTools.hxx>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

#include "conveniences.h"
#include "constants.h"

int main(int argc, char* argv[])
{

  std::string dirString = "/Users/ashraymalhotra/Desktop/Academic/VTKDev/Data/AllImages/";
  if (argc < 2)
    {
    cerr << "WARNING: Expected arguments" << argv[0] << "/directory/path/to/numbered/image/files/" << endl;
    cerr << "defaulting to " << dirString << endl;
    }
  else
    {
    dirString = argv[1];
    }

  // Have to define these variables as global since used in mapper, can be optimised better
  vtkImagePyramid *Pyramid[3];

  vtkImagePyramid *lowPass1[3];
  vtkImagePyramid *lowPass2[3];
  vtkImagePyramid *Difference[3];
  vtkImageData * FrameDifference[3];
  vtkImageData * OutputFrame[3];
  for (int i = 0; i < 3; i++)
    {
    FrameDifference[i] = vtkImageData::New();
    OutputFrame[i] = vtkImageData::New();
    }

  vtkSmartPointer<vtkImageExtractComponents> imageSource; // Moved outside since we need it to add back the differences
  vtkSmartPointer<vtkImageYIQToRGB> rgbConversionFilter =
    vtkSmartPointer<vtkImageYIQToRGB>::New(); // Made a global variable so that can be accessed by mapper



  int frameSize[NumberOfPyramidLevels];
  // Create reader to read all images
  vtksys::SystemTools::ConvertToUnixSlashes(dirString);
  vtkSmartPointer<vtkGlobFileNames> glob = vtkSmartPointer<vtkGlobFileNames>::New();
  glob->SetDirectory(dirString.c_str());
  glob->AddFileNames("*");
  int frameCount = glob->GetNumberOfFileNames();
  cout << "Number of Images read " << frameCount << "\n";

  for (int ImageNumber = 0; ImageNumber < frameCount; ImageNumber++)
    {
    std::string inputFilename = glob->GetNthFileName(ImageNumber);
    cout << inputFilename << "\n";

    // Read the file
    vtkSmartPointer<vtkImageReader2Factory> readerFactory =
      vtkSmartPointer<vtkImageReader2Factory>::New();
    vtkImageReader2 * imageReader = readerFactory->CreateImageReader2(inputFilename.c_str());
    imageReader->SetFileName(inputFilename.c_str());
    imageReader->Update();
    int* imageDimensionArray = imageReader->GetOutput()->GetExtent();
    int imageDimension1 = imageDimensionArray[1] - imageDimensionArray[0] + 1;
    int imageDimension2 = imageDimensionArray[3] - imageDimensionArray[2] + 1;

    // ---------------------Get YIQ components-----------------------------
    vtkSmartPointer<vtkImageRGBToYIQ> yiqFilter =
      vtkSmartPointer<vtkImageRGBToYIQ>::New();
    yiqFilter->SetInputConnection(imageReader->GetOutputPort());
    yiqFilter->Update();


    vtkSmartPointer<vtkImageData> YIQs[3];
    vtkSmartPointer<vtkImageExtractComponents> extractFilter =
      vtkSmartPointer<vtkImageExtractComponents>::New();
    for (int color_channel=0; color_channel<3; color_channel++)
      {
      extractFilter->SetInputConnection(yiqFilter->GetOutputPort());
      extractFilter->SetComponents(color_channel);
      extractFilter->Update();
      YIQs[color_channel] = extractFilter->GetOutput();
      }

    //-------------------------------------------------------------------------

    // ------------------Create image pyramids------------------------------
    for (int color_channel=0; color_channel<3; color_channel++)
      {
      Pyramid[color_channel] = new vtkImagePyramid(YIQs[color_channel], NumberOfPyramidLevels);
      }
    //-------------------------------------------------------------------------

    for (int color_channel=0; color_channel<3; color_channel++)
      {
      // -------------initialising lowPass with first frame--------------------------
      // Initialising the Previous frame pyramid for frame 1. We initialise it to the
      // first frame(which means first value of frame difference is going to be zero)
      if (ImageNumber==0)
        {
        lowPass1[color_channel] = new vtkImagePyramid();
        lowPass2[color_channel] = new vtkImagePyramid();
        lowPass1[color_channel]->ShallowCopy(Pyramid[color_channel]);
        lowPass2[color_channel]->ShallowCopy(Pyramid[color_channel]);
        }
      else
        {
        //---------Temporal IIR(Infinite impulse response) filtering-----------
        // -------Updating lowpass variable------
        vtkSmartPointer<vtkImageWeightedSum> sumFilter1 = vtkSmartPointer<vtkImageWeightedSum>::New();
        vtkSmartPointer<vtkImageWeightedSum> sumFilter2 = vtkSmartPointer<vtkImageWeightedSum>::New();
        sumFilter1->SetWeight(0,r1);
        sumFilter1->SetWeight(1,(1-r1));
        sumFilter2->SetWeight(0,r2);
        sumFilter2->SetWeight(1,(1-r2));
        for (int k=0; k<NumberOfPyramidLevels; k++)
          {
          sumFilter1->SetInputData(Pyramid[color_channel]->vtkImagePyramidData[k]);
          sumFilter1->AddInputData(lowPass1[color_channel]->vtkImagePyramidData[k]);
          sumFilter1->Update();
          lowPass1[color_channel]->vtkImagePyramidData[k]->ShallowCopy(sumFilter1->GetOutput());

          sumFilter2->SetInputData(Pyramid[color_channel]->vtkImagePyramidData[k]);
          sumFilter2->AddInputData(lowPass2[color_channel]->vtkImagePyramidData[k]);
          sumFilter2->Update();
          lowPass2[color_channel]->vtkImagePyramidData[k]->ShallowCopy(sumFilter2->GetOutput());
          }

        //------Updated lowpass variable-------

        //----Image Pyramid difference for IIR filtering-----
        vtkSmartPointer<vtkImageDifference> differenceFilter = vtkSmartPointer<vtkImageDifference>::New();
        differenceFilter->AllowShiftOff();
        differenceFilter->AveragingOff();
        differenceFilter->SetAllowShift(0);
        differenceFilter->SetThreshold(0);
        vtkSmartPointer<vtkImageMathematics> imageMath = vtkSmartPointer<vtkImageMathematics>::New();
        Difference[color_channel] = new vtkImagePyramid(NumberOfPyramidLevels);
        imageMath->SetOperationToSubtract();
        for (int k=0; k<NumberOfPyramidLevels; k++)
          {
          imageMath->SetInput1Data(lowPass1[color_channel]->vtkImagePyramidData[k]);
          imageMath->SetInput2Data(lowPass2[color_channel]->vtkImagePyramidData[k]);
          imageMath->Update();
          Difference[color_channel]->vtkImagePyramidData[k]->ShallowCopy(imageMath->GetOutput());
          }
        // ------End of pyramid difference------
        // ------------------End of temporal filtering---------------------------

        // ----------------Get image dimensions for spatial filtering------------
        for (int k=0; k<NumberOfPyramidLevels; k++)
          {
          int* a = lowPass1[color_channel]->vtkImagePyramidData[k]->GetExtent();
          int imageDimension1 = a[1] - a[0] + 1;
          int imageDimension2 = a[3] - a[2] + 1;
          frameSize[k] = pow((pow(imageDimension1,2) + pow(imageDimension2,2)), 0.5)/3;
          // 3 is an experimental constant used by the authors of the paper
          }
        // ----------------------------------------------------------------------

        vtkSmartPointer<vtkImageMathematics> differenceBooster = vtkSmartPointer<vtkImageMathematics>::New();
        differenceBooster->SetOperationToMultiplyByK();

        for (int k=1; k<NumberOfPyramidLevels-1; k++)
          {
          int currAlpha = frameSize[k]/(delta*8) - 1;
          currAlpha = currAlpha * exaggeration_factor;
          int mutiplier = currAlpha;
          if (mutiplier > alpha)
            {
            mutiplier = alpha;
            }
          // ----------Verify that multiplier is a float(or acceptable data format for SetConstantK)--------
          differenceBooster->SetConstantK(mutiplier);
          differenceBooster->SetInput1Data(Difference[color_channel]->vtkImagePyramidData[k]);
          differenceBooster->Update();
          Difference[color_channel]->vtkImagePyramidData[k]->ShallowCopy(differenceBooster->GetOutput());
          }
          // -------------End of spatial filtering----------------------

        //-----------Collapse the image Pyramid------------------------
        vtkSmartPointer<vtkImageResize> resize;
        resize = vtkSmartPointer<vtkImageResize>::New();
        resize->SetResizeMethodToOutputDimensions();

        vtkSmartPointer<vtkImageWeightedSum> sumFilter = vtkSmartPointer<vtkImageWeightedSum>::New();
        sumFilter->SetWeight(0,0.5);
        sumFilter->SetWeight(1,0.5);

        for (int g = (NumberOfPyramidLevels-1); g>=1; g--)
          {
          if (g == (NumberOfPyramidLevels-1))
            {
            resize->SetInputData(Difference[color_channel]->vtkImagePyramidData[g]);
            }
          else
            {
            resize->SetInputData(sumFilter->GetOutput());
            }

          resize->SetOutputDimensions(imageDimension1/pow(2,(g-1)), imageDimension2/pow(2,(g-1)), -1);
          resize->Update();

          vtkSmartPointer<vtkImageConvolve> convolveFilter2 =
            vtkSmartPointer<vtkImageConvolve>::New();
          convolveFilter2->SetInputData(resize->GetOutput());
          double kernel[25] = {1,4,6,4,1,
                               4,16,24,16,4,
                               6,24,36,24,6,
                               4,16,24,16,4,
                               1,4,6,4,1};
          for (int internal_loop=0; internal_loop<25;internal_loop++)
            {
            kernel[internal_loop] = kernel[internal_loop]/256;
            }
          convolveFilter2->SetKernel5x5(kernel);
          convolveFilter2->Update();
          sumFilter->SetInputData(convolveFilter2->GetOutput());
          sumFilter->AddInputData(Difference[color_channel]->vtkImagePyramidData[g-1]);
          sumFilter->Update();

          // --------------Save the final image in corresponding difference variable---------
          if (g==1)
            {
            FrameDifference[color_channel]->ShallowCopy(sumFilter->GetOutput());
            }
          }

        //---------------Pyramid Collapsed into image---------------------------

        //----------Chromatic Abberation to reduce noise---------------
        vtkSmartPointer<vtkImageMathematics> chromaticCorrection =
          vtkSmartPointer<vtkImageMathematics>::New();
        chromaticCorrection->SetOperationToMultiplyByK();
        chromaticCorrection->SetConstantK(chromatic_abberation);
        if (color_channel != YChannel)
          {
          chromaticCorrection->SetInput1Data(FrameDifference[color_channel]);
          chromaticCorrection->Update();
          FrameDifference[color_channel]->ShallowCopy(chromaticCorrection->GetOutput());
          }

        //--------Add back frame difference to the original frame that we have read------
        vtkSmartPointer<vtkImageWeightedSum> addDifferenceOrigFrameFilter =
          vtkSmartPointer<vtkImageWeightedSum>::New();
        // Note that we might have to multiply the intensity with a factor of 2 later...
        addDifferenceOrigFrameFilter->SetWeight(0,.5);
        addDifferenceOrigFrameFilter->SetWeight(1,.5);
        addDifferenceOrigFrameFilter->SetInputData(YIQs[color_channel]);
        addDifferenceOrigFrameFilter->AddInputData(FrameDifference[color_channel]);
        addDifferenceOrigFrameFilter->Update();
        OutputFrame[color_channel]->ShallowCopy(addDifferenceOrigFrameFilter->GetOutput());

        vtkSmartPointer<vtkImageMathematics> IntensityNormalisation =
          vtkSmartPointer<vtkImageMathematics>::New();
        IntensityNormalisation->SetOperationToMultiplyByK();
        IntensityNormalisation->SetConstantK(2);
        IntensityNormalisation->SetInput1Data(OutputFrame[color_channel]);
        IntensityNormalisation->Update();
        OutputFrame[color_channel]->ShallowCopy(IntensityNormalisation->GetOutput());
        //---------------------------------------------------------------------

        }   // End of the else loop(to perform operations only for frame numbers agreater than 1)
      } //End of iteration over the 3 color channels

    // ---------------------Combine color channels in 1 image------------------------

    if (ImageNumber!=0)
      {
      vtkSmartPointer<vtkImageAppendComponents> appendFilter =
        vtkSmartPointer<vtkImageAppendComponents>::New();
      //cerr << "Y dims: " << showDims(OutputFrame[YChannel]) << endl
      //     << "I dims: " << showDims(OutputFrame[IChannel]) << endl
      //     << "Q dims: " << showDims(OutputFrame[QChannel]) << endl;
      appendFilter->AddInputData(OutputFrame[YChannel]);
      appendFilter->AddInputData(OutputFrame[IChannel]);
      appendFilter->AddInputData(OutputFrame[QChannel]);
      appendFilter->Update();
      // -------------------------------------------------------------------------------

      // ---------------------Convert the YIQ frame to the RGB frame--------------------
      rgbConversionFilter->SetInputConnection(appendFilter->GetOutputPort());
      rgbConversionFilter->Update();
      // -------------------------------------------------------------------------------

      std::string iterationNumberString = to_string(ImageNumber);
      std::string outputFileName = "OutputFrame" + iterationNumberString+".png";
      vtkSmartPointer<vtkPNGWriter> writeDifferenceFrames = vtkSmartPointer<vtkPNGWriter>::New();
      writeDifferenceFrames->SetFileName(outputFileName.c_str());
      writeDifferenceFrames->SetInputData(rgbConversionFilter->GetOutput());
      writeDifferenceFrames->Write();
    }

  } //End of iteration over all input frames of the video(input images)
  return EXIT_SUCCESS;
}
