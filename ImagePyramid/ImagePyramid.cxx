#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkImageActor.h>
#include <vtkImageViewer2.h>
#include <vtkImageReader2Factory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageGaussianSmooth.h>

#include <vtkSmartPointer.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageReader2.h>
#include <vtkImageData.h>
#include <vtkImageMapper3D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkImageExtractComponents.h>


#include <vtkImageResize.h>
#include "vtkImageSincInterpolator.h"
#include <vtkImageRGBToYIQ.h>
#include <vtkImageYIQToRGB.h>

#include <vtkGlobFileNames.h>
#include <vtksys/SystemTools.hxx>
#include <string>
#include <sstream>
#include <cstdio>
#include <cmath>
#include "vtkStringArray.h"

#include <vtkImageShiftScale.h>
#include <vtkImageDifference.h>
#include <vtkImageMathematics.h>
#include <vtkImageWeightedSum.h>
#include <vtkImageAppendComponents.h>
#include <vtkPNGWriter.h>

std::string showDims (vtkImageData *img)
{
  std::ostringstream strm;
  int *exts;
  exts = img->GetExtent();
  strm << "(" << exts[1] - exts[0] << ", " << exts[3] - exts[2] << ", " << exts[5] - exts[4] << ")";
  return strm.str();
}

// Originally defined as struct but with a struct the variables go out of scope and hence get deleted.
class imageStruct
{
  public:
  imageStruct() { this->imagedata = vtkSmartPointer<vtkImageData>::New(); }
  ~imageStruct() { this->imagedata = NULL; }
  vtkSmartPointer<vtkImageData> imagedata;
};

int main(int argc, char* argv[])
{
  // Have to define these variables as global since used in mapper
  int NumberOfPyramidLevels = 6;
  imageStruct *imagePyramid = new imageStruct[NumberOfPyramidLevels];
  imageStruct *YPyramid = new imageStruct[NumberOfPyramidLevels];
  imageStruct *IPyramid = new imageStruct[NumberOfPyramidLevels];
  imageStruct *QPyramid = new imageStruct[NumberOfPyramidLevels];

  imageStruct *YPyramidPreviousFrame = new imageStruct[NumberOfPyramidLevels];
  imageStruct *IPyramidPreviousFrame = new imageStruct[NumberOfPyramidLevels];
  imageStruct *QPyramidPreviousFrame = new imageStruct[NumberOfPyramidLevels];
  imageStruct *YDifference = new imageStruct[NumberOfPyramidLevels];
  imageStruct *IDifference = new imageStruct[NumberOfPyramidLevels];
  imageStruct *QDifference = new imageStruct[NumberOfPyramidLevels];

  imageStruct *lowPass1Y = new imageStruct[NumberOfPyramidLevels];
  imageStruct *lowPass2Y = new imageStruct[NumberOfPyramidLevels];
  imageStruct *lowPass1I = new imageStruct[NumberOfPyramidLevels];
  imageStruct *lowPass2I = new imageStruct[NumberOfPyramidLevels];
  imageStruct *lowPass1Q = new imageStruct[NumberOfPyramidLevels];
  imageStruct *lowPass2Q = new imageStruct[NumberOfPyramidLevels];

  vtkSmartPointer<vtkImageData> FrameDifferenceY = vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkImageData> FrameDifferenceI = vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkImageData> FrameDifferenceQ = vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkImageData> OutputFrameY = vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkImageData> OutputFrameI = vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkImageData> OutputFrameQ = vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkImageExtractComponents> imageSource; // Moved outside since we need it to add back the differences
  vtkSmartPointer<vtkImageYIQToRGB> rgbConversionFilter =
    vtkSmartPointer<vtkImageYIQToRGB>::New(); // Made a global variable so that can be accessed by mapper



  int frameSize[NumberOfPyramidLevels];

  //TODO: Parse command line arguments
  std::string directoryName = "/Users/ashraymalhotra/Desktop/Academic/VTKDev/Data/ImageSequence/";

  // Create reader to read all images
  std::string dirString = directoryName;
  vtksys::SystemTools::ConvertToUnixSlashes(dirString);
  vtkSmartPointer<vtkGlobFileNames> glob = vtkSmartPointer<vtkGlobFileNames>::New();
  glob->SetDirectory(dirString.c_str());
  glob->AddFileNames("*");
  int frameCount = glob->GetNumberOfFileNames();
  cout << "Number of Images read " << frameCount << "\n";

  for (int ImageNumber = 0; ImageNumber < frameCount; ImageNumber++){
    std::string inputFilename = glob->GetNthFileName(ImageNumber);
    cout << inputFilename << "\n";

    // std::string inputFilename = "image.png";

    // Read the file
    vtkSmartPointer<vtkImageReader2Factory> readerFactory =
      vtkSmartPointer<vtkImageReader2Factory>::New();
    vtkImageReader2 * imageReader = readerFactory->CreateImageReader2(inputFilename.c_str());
    imageReader->SetFileName(inputFilename.c_str());
    imageReader->Update();
    // cout << __LINE__<<" Reached here \n";

    //--------------This part should go into an extract YIQ filter----------------
    vtkSmartPointer<vtkImageRGBToYIQ> yiqFilter =
      vtkSmartPointer<vtkImageRGBToYIQ>::New();
    yiqFilter->SetInputConnection(imageReader->GetOutputPort());
    yiqFilter->Update();

    vtkSmartPointer<vtkImageExtractComponents> extractYFilter =
      vtkSmartPointer<vtkImageExtractComponents>::New();
    extractYFilter->SetInputConnection(yiqFilter->GetOutputPort());
    extractYFilter->SetComponents(0);
    extractYFilter->Update();

    vtkSmartPointer<vtkImageExtractComponents> extractIFilter =
      vtkSmartPointer<vtkImageExtractComponents>::New();
    extractIFilter->SetInputConnection(yiqFilter->GetOutputPort());
    extractIFilter->SetComponents(1);
    extractIFilter->Update();

    vtkSmartPointer<vtkImageExtractComponents> extractQFilter =
      vtkSmartPointer<vtkImageExtractComponents>::New();
    extractQFilter->SetInputConnection(yiqFilter->GetOutputPort());
    extractQFilter->SetComponents(2);
    extractQFilter->Update();
    //-----------------------------------------------------------------------------

    int NumberOfPyramidLevels = 6;
    for (int j=0; j<3; j++){
      switch (j){
        case 0:
          imageSource = extractYFilter;
          break;
        case 1:
          imageSource = extractIFilter;
          break;
        case 2:
          imageSource = extractQFilter;
          break;
      }

      // ---------------------Image Pyramid Construction---------------------
      int* a = imageReader->GetOutput()->GetExtent();
      int imageDimension1 = a[1] - a[0] + 1;
      int imageDimension2 = a[3] - a[2] + 1;

      vtkSmartPointer<vtkImageGaussianSmooth> gaussianSmoothFilter;
      vtkSmartPointer<vtkImageResize> resize;

      imagePyramid[0].imagedata->ShallowCopy(imageSource->GetOutput());

      // if (!imagePyramid[0].imagedata) printf("0 is NULL!\n");
      for (int i=1; i<NumberOfPyramidLevels; i++){
        gaussianSmoothFilter = vtkSmartPointer<vtkImageGaussianSmooth>::New();
        // gaussianSmoothFilter->SetInputData(imagePyramid[i-1].imagedata);
        gaussianSmoothFilter->SetInputData(imagePyramid[i-1].imagedata);
        gaussianSmoothFilter->Update();

        // Can set standard deviation of the gaussian filter using "SetStandardDeviation (double std)"

        resize = vtkSmartPointer<vtkImageResize>::New();
        resize->SetResizeMethodToOutputDimensions();
#if VTK_MAJOR_VERSION <= 5
        resize->SetInput(gaussianSmoothFilter->GetOutput());
#else
        resize->SetInputData(gaussianSmoothFilter->GetOutput());
#endif
        resize->SetOutputDimensions(imageDimension1/(pow(2,i)), imageDimension2/(pow(2,i)), 1);
        resize->Update();
        imagePyramid[i].imagedata->ShallowCopy(resize->GetOutput());
      }
      // ------------------Image Pyramid construction complete--------------

      // -------------Copy Image Pyramid into corresponding variable-------
      switch (j){
        case 0:
          {
            //Make the code below in a seperate function to copy image pyramids
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              YPyramid[k].imagedata->ShallowCopy(imagePyramid[k].imagedata);
            }
            break;
          }
        case 1:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              IPyramid[k].imagedata->ShallowCopy(imagePyramid[k].imagedata);
            }
            break;
          }
        case 2:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              QPyramid[k].imagedata->ShallowCopy(imagePyramid[k].imagedata);
            }
            break;
          }
      }
      // --------Copied image pyramid into corresponding variable-------------

      float r1 = 0.4;
      float r2 = 0.05;

      // -------------initialising lowPass with first frame--------------------------
      // Initialising the Previous frame pyramid for frame 1. We initialise it to the first frame(which means first value of frame difference is going to be zero)
      if (ImageNumber==0)
      {
        switch (j)
        {
          case 0:
          {
            //Make the code below in a seperate function to copy image pyramids
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              lowPass1Y[k].imagedata->ShallowCopy(YPyramid[k].imagedata);
              lowPass2Y[k].imagedata->ShallowCopy(YPyramid[k].imagedata);
            }
            break;
          }
          case 1:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              lowPass1I[k].imagedata->ShallowCopy(IPyramid[k].imagedata);
              lowPass2I[k].imagedata->ShallowCopy(IPyramid[k].imagedata);
            }
            break;
          }
          case 2:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              lowPass1Q[k].imagedata->ShallowCopy(QPyramid[k].imagedata);
              lowPass2Q[k].imagedata->ShallowCopy(QPyramid[k].imagedata);
            }
            break;
          }
        }
      } //--------------------------------------------------------------------
      else
      {
        //---------Temporal IIR(Infinite impulse response) filtering-----------
        // ----Updating lowpass variable------
        vtkSmartPointer<vtkImageWeightedSum> sumFilter1 = vtkSmartPointer<vtkImageWeightedSum>::New();
        vtkSmartPointer<vtkImageWeightedSum> sumFilter2 = vtkSmartPointer<vtkImageWeightedSum>::New();
        sumFilter1->SetWeight(0,r1);
        sumFilter1->SetWeight(1,(1-r1));
        sumFilter2->SetWeight(0,r2);
        sumFilter2->SetWeight(1,(1-r2));
        switch (j)
        {
          case 0:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              sumFilter1->SetInputData(YPyramid[k].imagedata);
              sumFilter1->AddInputData(lowPass1Y[k].imagedata);
              sumFilter1->Update();
              lowPass1Y[k].imagedata->ShallowCopy(sumFilter1->GetOutput());
              sumFilter2->SetInputData(YPyramid[k].imagedata);
              sumFilter2->AddInputData(lowPass2Y[k].imagedata);
              sumFilter2->Update();
              lowPass2Y[k].imagedata->ShallowCopy(sumFilter2->GetOutput());
            }
            break;
          }
          case 1:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              sumFilter1->SetInputData(IPyramid[k].imagedata);
              sumFilter1->AddInputData(lowPass1I[k].imagedata);
              sumFilter1->Update();
              lowPass1I[k].imagedata->ShallowCopy(sumFilter1->GetOutput());
              sumFilter2->SetInputData(IPyramid[k].imagedata);
              sumFilter2->AddInputData(lowPass2I[k].imagedata);
              sumFilter2->Update();
              lowPass2I[k].imagedata->ShallowCopy(sumFilter2->GetOutput());
            }
            break;
          }
          case 2:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              sumFilter1->SetInputData(QPyramid[k].imagedata);
              sumFilter1->AddInputData(lowPass1Q[k].imagedata);
              sumFilter1->Update();
              lowPass1Q[k].imagedata->ShallowCopy(sumFilter1->GetOutput());
              sumFilter2->SetInputData(QPyramid[k].imagedata);
              sumFilter2->AddInputData(lowPass2Q[k].imagedata);
              sumFilter2->Update();
              lowPass2Q[k].imagedata->ShallowCopy(sumFilter2->GetOutput());
            }
            break;
          }
        }
        //-----Updated lowpass variable-------

        //----Image Pyramid difference for IIR filtering------
        vtkSmartPointer<vtkImageDifference> differenceFilter = vtkSmartPointer<vtkImageDifference>::New();
        differenceFilter->AllowShiftOff();
        differenceFilter->AveragingOff();
        differenceFilter->SetAllowShift(0);
        differenceFilter->SetThreshold(0);
        // TODO Verify that image mathematics works the way we would expect it to!!
        vtkSmartPointer<vtkImageMathematics> imageMath = vtkSmartPointer<vtkImageMathematics>::New();
        imageMath->SetOperationToSubtract();
        switch(j)
        {
          case 0:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              imageMath->SetInput1Data(lowPass1Y[k].imagedata);
              imageMath->SetInput2Data(lowPass2Y[k].imagedata);
              imageMath->Update();
              YDifference[k].imagedata->ShallowCopy(imageMath->GetOutput());
            }
            break;
          }
          case 1:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              imageMath->SetInput1Data(lowPass1I[k].imagedata);
              imageMath->SetInput2Data(lowPass2I[k].imagedata);
              imageMath->Update();
              IDifference[k].imagedata->ShallowCopy(imageMath->GetOutput());
            }
            break;
          }
          case 2:
          {
            for (int k=0; k<NumberOfPyramidLevels; k++)
            {
              imageMath->SetInput1Data(lowPass1Q[k].imagedata);
              imageMath->SetInput2Data(lowPass2Q[k].imagedata);
              imageMath->Update();
              QDifference[k].imagedata->ShallowCopy(imageMath->GetOutput());
            }
            break;
          }
        }
        // ------End of pyramid difference------
      // ------------------End of temporal filtering-------------------------------------------

      // ----------------Get image dimensions for spatial filtering----------------------
        for (int k=0; k<NumberOfPyramidLevels; k++)
        {
          // Will the below code work or not?? Is imagedata returned as a pointer? Verify.
          int* a = lowPass1Q[k].imagedata->GetExtent();
          int imageDimension1 = a[1] - a[0] + 1;
          int imageDimension2 = a[3] - a[2] + 1;
          frameSize[k] = pow((pow(imageDimension1,2) + pow(imageDimension2,2)), 0.5)/3;
          // 3 is an experimental constant used by the authors of the paper
        }
      // ---------------------------------------------------------------------------------


      // ----------------Spatial Filtering----------------------
      // Refer to the diagram included in the paper for explanation of the code below

      // ----------Code below causes seg fault--------------------

      //---------Assign the top and bottom of the image pyramid to zero--------
        // switch(j)
        // {
        //   case 0:
        //   {
        //     int* a = YDifference[0].imagedata->GetExtent();
        //     int Y0Dimension1 = a[1] - a[0] + 1;
        //     int Y0Dimension2 = a[3] - a[2] + 1;
        //     for (int y = 0; y < Y0Dimension2; y++)
        //     {
        //       for (int x = 0; x < Y0Dimension1; x++)
        //       {
        //       double* pixel = static_cast<double*>(YDifference[0].imagedata->GetScalarPointer(x,y,0));
        //       pixel[0] = 0.0;
        //       }
        //     }

        //     a = YDifference[NumberOfPyramidLevels-1].imagedata->GetExtent();
        //     int YNDimension1 = a[1] - a[0] + 1;
        //     int YNDimension2 = a[3] - a[2] + 1;
        //     for (int y = 0; y < YNDimension2; y++)
        //     {
        //       for (int x = 0; x < YNDimension1; x++)
        //       {
        //       double* pixel = static_cast<double*>(YDifference[NumberOfPyramidLevels-1].imagedata->GetScalarPointer(x,y,0));
        //       pixel[0] = 0.0;
        //       }
        //     }
        //   break;
        //   }
        //   case 1:
        //   {
        //     int* a = IDifference[0].imagedata->GetExtent();
        //     int Y0Dimension1 = a[1] - a[0] + 1;
        //     int Y0Dimension2 = a[3] - a[2] + 1;
        //     for (int y = 0; y < Y0Dimension2; y++)
        //     {
        //       for (int x = 0; x < Y0Dimension1; x++)
        //       {
        //       double* pixel = static_cast<double*>(IDifference[0].imagedata->GetScalarPointer(x,y,0));
        //       pixel[0] = 0.0;
        //       }
        //     }

        //     a = IDifference[NumberOfPyramidLevels-1].imagedata->GetExtent();
        //     int YNDimension1 = a[1] - a[0] + 1;
        //     int YNDimension2 = a[3] - a[2] + 1;
        //     for (int y = 0; y < YNDimension2; y++)
        //     {
        //       for (int x = 0; x < YNDimension1; x++)
        //       {
        //       double* pixel = static_cast<double*>(IDifference[NumberOfPyramidLevels-1].imagedata->GetScalarPointer(x,y,0));
        //       pixel[0] = 0.0;
        //       }
        //     }
        //   break;
        //   }
        //   case 2:
        //   {
        //     int* a = QDifference[0].imagedata->GetExtent();
        //     int Y0Dimension1 = a[1] - a[0] + 1;
        //     int Y0Dimension2 = a[3] - a[2] + 1;
        //     for (int y = 0; y < Y0Dimension2; y++)
        //     {
        //       for (int x = 0; x < Y0Dimension1; x++)
        //       {
        //       double* pixel = static_cast<double*>(QDifference[0].imagedata->GetScalarPointer(x,y,0));
        //       pixel[0] = 0.0;
        //       }
        //     }

        //     a = QDifference[NumberOfPyramidLevels-1].imagedata->GetExtent();
        //     int YNDimension1 = a[1] - a[0] + 1;
        //     int YNDimension2 = a[3] - a[2] + 1;
        //     for (int y = 0; y < YNDimension2; y++)
        //     {
        //       for (int x = 0; x < YNDimension1; x++)
        //       {
        //       double* pixel = static_cast<double*>(QDifference[NumberOfPyramidLevels-1].imagedata->GetScalarPointer(x,y,0));
        //       pixel[0] = 0.0;
        //       }
        //     }
        //   break;
        //   }
        // }
        // --------------------------------------------------------------

        int alpha = 10;
        int lambda_c = 16;
        float delta;
        delta = lambda_c/8/(1+alpha);
        //exaggeration_factor HAS to be a user defined constant, above values could be hardcoded as well
        int exaggeration_factor = 20;
        vtkSmartPointer<vtkImageMathematics> differenceBooster = vtkSmartPointer<vtkImageMathematics>::New();
        differenceBooster->SetOperationToMultiplyByK();

        switch(j)
        {
          case 0:
          {
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
              // differenceBooster->SetInputConnection(imageSource->GetOutputPort());
              differenceBooster->SetInput1Data(YDifference[k].imagedata);
              differenceBooster->Update();
              YDifference[k].imagedata->ShallowCopy(differenceBooster->GetOutput());
            }
            break;
          }
          case 1:
          {
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
              // differenceBooster->SetInputConnection(imageSource->GetOutputPort());
              differenceBooster->SetInput1Data(IDifference[k].imagedata);
              differenceBooster->Update();
              IDifference[k].imagedata->ShallowCopy(differenceBooster->GetOutput());
            }
            break;
          }
          case 2:
          {
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
              // differenceBooster->SetInputConnection(imageSource->GetOutputPort());
              differenceBooster->SetInput1Data(QDifference[k].imagedata);
              differenceBooster->Update();
              QDifference[k].imagedata->ShallowCopy(differenceBooster->GetOutput());
            }
            break;
          }
        }
        // -------------End of spatial filtering----------------------

        //-------------VERIFY THE ALGO TO COLLAPSE PYRAMID. SEEMS WRONG---------------
        //-----------Collapse the image Pyramid------------------------

        // Verify that the algorithm to implement the laplacian pyramid is correct
        resize = vtkSmartPointer<vtkImageResize>::New();
        gaussianSmoothFilter = vtkSmartPointer<vtkImageGaussianSmooth>::New();
        vtkSmartPointer<vtkImageWeightedSum> sumFilter = vtkSmartPointer<vtkImageWeightedSum>::New();
        sumFilter->SetWeight(0,0.5);
        sumFilter->SetWeight(1,0.5);
        resize->SetResizeMethodToOutputDimensions();
        for (int g = (NumberOfPyramidLevels-1); g>=1; g--)
        {
          if (g == (NumberOfPyramidLevels-1))
          {
            switch (j) {
              case 0:
              {
                #if VTK_MAJOR_VERSION <= 5
                          resize->SetInput(YDifference[g].imagedata);
                #else
                          resize->SetInputData(YDifference[g].imagedata);
                #endif
                break;
              }
              case 1:
              {
                #if VTK_MAJOR_VERSION <= 5
                          resize->SetInput(IDifference[g].imagedata);
                #else
                          resize->SetInputData(IDifference[g].imagedata);
                #endif
                break;
              }
              case 2:
              {
                #if VTK_MAJOR_VERSION <= 5
                          resize->SetInput(QDifference[g].imagedata);
                #else
                          resize->SetInputData(QDifference[g].imagedata);
                #endif
                break;
              }
            }
          }
          else
          {
  #if VTK_MAJOR_VERSION <= 5
            resize->SetInput(sumFilter->GetOutput());
  #else
            resize->SetInputData(sumFilter->GetOutput());
  #endif
          }

          sumFilter = NULL;
          sumFilter = vtkSmartPointer<vtkImageWeightedSum>::New();
          sumFilter->SetWeight(0,0.5);
          sumFilter->SetWeight(1,0.5);

          resize->SetOutputDimensions(imageDimension1/pow(2,(g-1)), imageDimension2/pow(2,(g-1)), -1);
          resize->Update();

          gaussianSmoothFilter->SetInputData(resize->GetOutput());
          gaussianSmoothFilter->Update();
          sumFilter->SetInputData(gaussianSmoothFilter->GetOutput());
          switch (j) {
            case 0:
            {
              sumFilter->AddInputData(YDifference[g-1].imagedata);
              break;
            }
            case 1:
            {
              sumFilter->AddInputData(IDifference[g-1].imagedata);
              break;
            }
            case 2:
            {
              sumFilter->AddInputData(QDifference[g-1].imagedata);
              break;
            }
          }
          sumFilter->Update();

          // --------------Save the final image in corresponding difference variable---------
          if (g==1)
          {
            switch (j)
            {
              case 0:
              {
                FrameDifferenceY->ShallowCopy(sumFilter->GetOutput());
              }
              case 1:
              {
                FrameDifferenceI->ShallowCopy(sumFilter->GetOutput());
                break;
              }
              case 2:
              {
                FrameDifferenceQ->ShallowCopy(sumFilter->GetOutput());
                break;
              }
            }
          }
        }
        //---------------Pyramid Collapsed into image---------------------------



        //----------Chromatic Abberation to reduce noise---------------
        double chromatic_abberation = 0.1; //User given input
        vtkSmartPointer<vtkImageMathematics> chromaticCorrection =
          vtkSmartPointer<vtkImageMathematics>::New();
        chromaticCorrection->SetOperationToMultiplyByK();
        chromaticCorrection->SetConstantK(chromatic_abberation);
        switch(j)
        {
          // Do nothing for Y channel
          case 1:
          {
            chromaticCorrection->SetInput1Data(FrameDifferenceI);
            chromaticCorrection->Update();
            FrameDifferenceI->ShallowCopy(chromaticCorrection->GetOutput());
            break;
          }
          case 2:
          {
            chromaticCorrection->SetInput1Data(FrameDifferenceQ);
            chromaticCorrection->Update();
            FrameDifferenceQ->ShallowCopy(chromaticCorrection->GetOutput());
          }
        }

        //--------Add back frame difference to the original frame that we have read------
        vtkSmartPointer<vtkImageWeightedSum> addDifferenceOrigFrameFilter =
        vtkSmartPointer<vtkImageWeightedSum>::New();
        // Note that we might have to multiply the intensity with a factor of 2 later...
        addDifferenceOrigFrameFilter->SetWeight(0,.5);
        addDifferenceOrigFrameFilter->SetWeight(1,.5);
        // addDifferenceOrigFrameFilter->SetInputData(imageSource->GetOutput());
        switch (j) {
          case 0:
          {
            addDifferenceOrigFrameFilter->SetInputData(extractYFilter->GetOutput());
            addDifferenceOrigFrameFilter->AddInputData(FrameDifferenceY);
            addDifferenceOrigFrameFilter->Update();
            OutputFrameY->ShallowCopy(addDifferenceOrigFrameFilter->GetOutput());
            break;
          }
          case 1:
          {
            addDifferenceOrigFrameFilter->SetInputData(extractIFilter->GetOutput());
            addDifferenceOrigFrameFilter->AddInputData(FrameDifferenceI);
            addDifferenceOrigFrameFilter->Update();
            OutputFrameI->ShallowCopy(addDifferenceOrigFrameFilter->GetOutput());
            break;
          }
          case 2:
          {
            addDifferenceOrigFrameFilter->SetInputData(extractQFilter->GetOutput());
            addDifferenceOrigFrameFilter->AddInputData(FrameDifferenceQ);
            addDifferenceOrigFrameFilter->Update();
            OutputFrameQ->ShallowCopy(addDifferenceOrigFrameFilter->GetOutput());
            break;
          }
        }
        //---------------------------------------------------------------------

        // -----------------Debugging tip. Try to look at the color channel differences. We know what to expect in each of the color channels

      }   // End of the else loop(to perform operations only for frame numbers greater than 1)
    } //End of iteration over the 3 color channels

    // ---------------------Combine color channels in 1 image------------------------

    if(ImageNumber!=0)
    {
      vtkSmartPointer<vtkImageAppendComponents> appendFilter =
        vtkSmartPointer<vtkImageAppendComponents>::New();
      cerr << "Y dims: " << showDims(OutputFrameY) << endl
           << "I dims: " << showDims(OutputFrameI) << endl
           << "Q dims: " << showDims(OutputFrameQ) << endl;
      appendFilter->AddInputData(OutputFrameY);
      appendFilter->AddInputData(OutputFrameI);
      appendFilter->AddInputData(OutputFrameQ);
      appendFilter->Update();
      // -------------------------------------------------------------------------------

      // ---------------------Convert the YIQ frame to the RGB frame--------------------
      rgbConversionFilter->SetInputConnection(appendFilter->GetOutputPort());
      rgbConversionFilter->Update();
      // -------------------------------------------------------------------------------

      // TODO - Output that frame in the mapper(or write the frame as a png file)

      std::string iterationNumberString = std::to_string(ImageNumber);
      std::string outputFileName = "OutputFrame" + iterationNumberString+".png";
      vtkSmartPointer<vtkPNGWriter> writeDifferenceFrames = vtkSmartPointer<vtkPNGWriter>::New();
      writeDifferenceFrames->SetFileName(outputFileName.c_str());
      // writeDifferenceFrames->SetInputConnection(scaleDifference->GetOutputPort());
      writeDifferenceFrames->SetInputData(rgbConversionFilter->GetOutput());
      // writeDifferenceFrames->SetInputData(OutputFrameY);
      writeDifferenceFrames->Write();
    }
    

  // vtkSmartPointer<vtkDataSetMapper> mapper =
  //   vtkSmartPointer<vtkDataSetMapper>::New();
  //
  // // mapper->SetInputConnection(resize->GetOutputPort());
  //
  // //-------------------Suspected code snippet causing segmentation error----------------------------
  // // Why does mapper->SetInputData(imagePyramid[0].imagedata);  throw up an error?
  // // mapper->SetInputData(YDifference[0].imagedata);
  // mapper->SetInputData(rgbConversionFilter->GetOutput());
  //   // mapper->SetInputData(differenceFilter->GetOutput);
  // // mapper->SetInputData(resize->GetOutput());
  // //--------------------------------------------------------------------------------------------------
  //
  // vtkSmartPointer<vtkActor> actor =
  //   vtkSmartPointer<vtkActor>::New();
  // actor->SetMapper(mapper);
  // actor->GetProperty()->SetRepresentationToWireframe();
  //
  // vtkSmartPointer<vtkRenderer> renderer =
  //   vtkSmartPointer<vtkRenderer>::New();
  // renderer->AddActor(actor);
  // renderer->ResetCamera();
  // renderer->SetBackground(1,1,1);
  //
  // vtkSmartPointer<vtkRenderWindow> renderWindow =
  //   vtkSmartPointer<vtkRenderWindow>::New();
  // renderWindow->AddRenderer(renderer);
  //
  // vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  //   vtkSmartPointer<vtkRenderWindowInteractor>::New();
  // renderWindowInteractor->SetRenderWindow(renderWindow);
  // renderWindowInteractor->Initialize();
  //
  // renderWindowInteractor->Start();
  } //End of iteration over all input frames of the video(input images)
  return EXIT_SUCCESS;
}