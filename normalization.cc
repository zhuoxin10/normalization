#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif



// Software Guide : EndCodeSnippet
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkListSample.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkNumericSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkNeighborhoodIterator.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkNumericSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkNeighborhoodIterator.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkAffineTransform.h"
#include "itkCastImageFilter.h"
#include "itkImportImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkGaussianBlurImageFunction.h>
#include <itkStatisticsImageFilter.h>
#include <itkChangeInformationImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include "itkImageDuplicator.h"
#include "itkPNGImageIOFactory.h"


#include <cstdlib>

#include <algorithm>
#include <functional>
#include <iostream>

#include "math.h"
#include <vector>
#include <iterator>
#include <string>
#include "utils.h"
#include <stdio.h>
#include <iostream>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>


using namespace std;

typedef unsigned short    PixelType;
const unsigned int      Dimension = 2;
typedef itk::Image< PixelType, Dimension >         ImageType;
typedef ImageType Flima;
ImageType::IndexType idx2;

std::vector<float> const sigmaV4All {4.0f};


auto const blur = [](itk::ImageSource<Flima> const& in, unsigned int axis, double sigma){
    typedef itk::RecursiveGaussianImageFilter<Flima, Flima> RecursiveGaussian;
    using Blur = RecursiveGaussian;
    auto bluh = Blur::New();
    bluh->SetInput(in.GetOutput());
    bluh->SetDirection(axis);
    bluh->SetSigma(sigma);
    bluh->SetNormalizeAcrossScale(true);
    return bluh;
};


auto const subtract = [](itk::ImageSource<Flima> const& lh, itk::ImageSource<Flima> const& rh){
    typedef itk::SubtractImageFilter <Flima, Flima> Subtract;
    auto shub = Subtract::New();
    shub->SetInput1(lh.GetOutput());
    shub->SetInput2(rh.GetOutput());
    return shub;
};



Flima::Pointer combineEenergyBands2Image(std::vector< Flima::Pointer >  energyBandsSource,
                                                     Flima::Pointer &maskSource){
    typedef itk::ImageDuplicator< ImageType > DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(maskSource);
    duplicator->Update();
    ImageType::Pointer image = duplicator->GetOutput();


    ImageType::SizeType sizeSource = maskSource->GetLargestPossibleRegion().GetSize();
//    for(size_t k = 0; k< sizeSource[2]; k++) {
        for(size_t j = 0; j< sizeSource[1]; j++) {
            for(size_t i = 0; i< sizeSource[0]; i++)  {

                idx2[0]=i;
                idx2[1]=j;
//                idx3[2]=k;

                if (maskSource->GetPixel(idx2)>0){

                float valueAtIJ=0;

                for (unsigned int e=0; e<energyBandsSource.size();e++) {
                    valueAtIJ+=energyBandsSource[e].GetPointer()->GetPixel(idx2);
                }

                image->SetPixel(idx2,valueAtIJ);
            }
            else{
                image->SetPixel(idx2,0);
            }
            }
        }
//    }
    return image;
}



std::vector< Flima::Pointer >  getEnergyBands( ImageType::Pointer  &image, std::vector<float> sigmaV) {

    typedef itk::CastImageFilter< Flima,Flima> CastFloat2Float;
    auto inputFloat = CastFloat2Float::New();

    std::vector< Flima::Pointer > energyBandsImgs; // the first one is the original image

    Flima::Pointer lastLPBand=image;

    for (size_t i =0; i<sigmaV.size(); i++){
        float sigma=sigmaV[i];
        inputFloat->SetInput(lastLPBand);

        auto blurred_x  = blur(*inputFloat, 0, sigma);
        auto blurred_xy = blur(*blurred_x , 1, sigma);
//        auto blurred_xyz = blur(*blurred_xy , 2, sigma);
        auto highPass = subtract(*inputFloat, *blurred_xy);

        highPass->Update();

        lastLPBand=blurred_xy->GetOutput();

        energyBandsImgs.push_back(highPass->GetOutput());

        if (i==(sigmaV.size()-1)){
            energyBandsImgs.push_back(blurred_xy->GetOutput());
        }

    }

    return energyBandsImgs;
}


float getMean(Flima::Pointer energyBandImg, ImageType::Pointer &mask){

    float totalV=0;
    float numV=0;

    ImageType::SizeType size = mask->GetLargestPossibleRegion().GetSize();
//    for(size_t k = 0; k< size[2]; k++) {
        for(size_t j = 0; j< size[1]; j++) {
            for(size_t i = 0; i< size[0]; i++)  {

                idx2[0]=i;
                idx2[1]=j;
//                idx3[2]=k;

                if (mask->GetPixel(idx2)>0){

                totalV+=energyBandImg.GetPointer()->GetPixel(idx2);
                numV++;
                }
            }
        }
//    }
    return totalV/numV;
}

float getStd(Flima::Pointer energyBandImg, ImageType::Pointer &mask, float mean){

    float totalV=0;
    float numV=0;

    ImageType::SizeType size = mask->GetLargestPossibleRegion().GetSize();
//    for(size_t k = 0; k< size[2]; k++) {
        for(size_t j = 0; j< size[1]; j++) {
            for(size_t i = 0; i< size[0]; i++)  {

                idx2[0]=i;
                idx2[1]=j;
//                idx3[2]=k;

                if (mask->GetPixel(idx2)>0){
                    totalV+=pow(energyBandImg.GetPointer()->GetPixel(idx2)-mean,2.0);
                    numV++;
                }
            }
        }
//    }
    return sqrt(totalV/numV-1);
}


void modifyEnergyBands(std::vector< Flima::Pointer >  &energyBandsSource,
                       std::vector<float> energyBandImgSourceMean,
                       std::vector<float> energyBandImgSourceStd,
                       std::vector<float> energyBandImgTargetMean,
                       std::vector<float> energyBandImgTargetStd,
                       Flima::Pointer &maskSource){
    Flima::IndexType index;

    for (unsigned int e=0; e<energyBandsSource.size();e++) {
        float stdSource=energyBandImgSourceStd[e];
        float stdTarget=energyBandImgTargetStd[e];

        float meanSource=energyBandImgSourceMean[e];
        float meanTarget=energyBandImgTargetMean[e];

        ImageType::SizeType sizeSource = maskSource->GetLargestPossibleRegion().GetSize();
//        for(size_t k = 0; k< sizeSource[2]; k++) {
            for(size_t j = 0; j< sizeSource[1]; j++) {
                for(size_t i = 0; i< sizeSource[0]; i++)  {

                    idx2[0]=i;
                    idx2[1]=j;
//                    idx3[2]=k;

                    if (maskSource->GetPixel(idx2)>0){

                        float pixelOrg=energyBandsSource[e].GetPointer()->GetPixel(idx2);

                        if (!std::isnan(stdSource)){
                            float pixelNew = (pixelOrg-meanSource)*stdTarget/stdSource + meanTarget ;
                            energyBandsSource[e].GetPointer()->SetPixel(idx2,pixelNew);
                        }
                    }

                }
            }
//        }

        /*std::cout<<"band "<<e<<std::endl;

        std::cout<<"meanSource "<<meanSource<<std::endl;
        std::cout<<"stdSource "<<stdSource<<std::endl;

        std::cout<<"meanTarget "<<meanTarget<<std::endl;
        std::cout<<"stdTarget "<<stdTarget<<std::endl;*/

    }
}


void matchingTwoImages( const string outname,  ImageType::Pointer &imageSouce, ImageType::Pointer  &imageTarget){

    auto normImg = imageSouce;

    typedef itk::ImageDuplicator< ImageType > DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(imageSouce);

    duplicator->Update();
    ImageType::Pointer maskSource = duplicator->GetOutput();

    duplicator->SetInputImage(imageTarget);
    duplicator->Update();
    ImageType::Pointer maskTarget = duplicator->GetOutput();


    ImageType::SizeType sizeSource = imageSouce->GetLargestPossibleRegion().GetSize();

    std::cout<<"sizeSource haha "<<sizeSource<<std::endl;

//    for(size_t k = 0; k< sizeSource[2]; k++) {
        for(size_t j = 0; j< sizeSource[1]; j++) {
            for(size_t i = 0; i< sizeSource[0]; i++)  {

                idx2[0]=i;
                idx2[1]=j;
//                idx3[2]=k;

                maskSource->SetPixel(idx2,1);
                if (imageSouce->GetPixel(idx2)<600){
                    maskSource->SetPixel(idx2,0);
                }else{
                    maskSource->SetPixel(idx2,1);
                }
            }
        }
//    }


    ImageType::SizeType sizeTarget = imageTarget->GetLargestPossibleRegion().GetSize();
//    for(size_t k = 0; k< sizeTarget[2]; k++) {
        for(size_t j = 0; j< sizeTarget[1]; j++) {
            for(size_t i = 0; i< sizeTarget[0]; i++)  {

                idx2[0]=i;
                idx2[1]=j;
//                idx3[2]=k;

                maskTarget->SetPixel(idx2,1);
                if (imageTarget->GetPixel(idx2)<30){
                    maskTarget->SetPixel(idx2,0);
                }else{
                    maskTarget->SetPixel(idx2,1);
                }
            }
        }
//    }


    std::vector< Flima::Pointer >  energyBandsTarget=getEnergyBands(imageTarget,sigmaV4All);

    std::vector<float> energyBandImgTargetMean;
    std::vector<float> energyBandImgTargetStd;
    for (unsigned int i = 0; i<energyBandsTarget.size(); i++){
        float meanTarget=getMean(energyBandsTarget[i], maskTarget);
        float stdTarget=getStd(energyBandsTarget[i], maskTarget, meanTarget);
        energyBandImgTargetMean.push_back(meanTarget);
        energyBandImgTargetStd.push_back(stdTarget);
    }

    // 5 iterations
    for (int it = 0; it < 5; it++) {
        std::cout<<"iteration "<<it<<std::endl;
        std::vector< Flima::Pointer >  energyBandsSource=getEnergyBands(normImg,sigmaV4All);
        std::vector<float> energyBandImgSourceMean;
        std::vector<float> energyBandImgSourceStd;
        for (unsigned int i =0; i<energyBandsSource.size(); i++){
            float meanSource=getMean(energyBandsSource[i], maskSource);
            float stdSource=getStd(energyBandsSource[i], maskSource, meanSource);
            energyBandImgSourceMean.push_back(meanSource);
            energyBandImgSourceStd.push_back(stdSource);
        }

        modifyEnergyBands(energyBandsSource,
                          energyBandImgSourceMean,
                          energyBandImgSourceStd,
                          energyBandImgTargetMean,
                          energyBandImgTargetStd,
                          maskSource);

        normImg=combineEenergyBands2Image(energyBandsSource, maskSource);

    }

    std::cout << "Energy bands target are: " << std::endl;
    std::cout << "- Mean: [";
    for (unsigned int i =0; i<energyBandImgTargetMean.size(); i++){
        std::cout << " " << energyBandImgTargetMean[i] << " ,";
    }
    std::cout << " ]" << std::endl;
    std::cout << "- Standard Dev: [";
    for (unsigned int i =0; i<energyBandImgTargetStd.size(); i++){
        std::cout << " " << energyBandImgTargetStd[i] << " ,";
    }
    std::cout << " ]" << std::endl;

    using WriterType = itk::ImageFileWriter< ImageType  >;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outname); //normalization后的图片
    writer->SetInput(normImg);

    /*try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error: " << error << std::endl;
      }

    writer->SetFileName("Souce.png");
    writer->SetInput(imageSouce);

    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error: " << error << std::endl;
      }

    writer->SetFileName("target.png");
    writer->SetInput(imageTarget);
    */

    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error: " << error << std::endl;
      }
}




//int readITKImage( const string filename, ImageType::Pointer &image ){
//
//
//
//    typedef itk::ImageSeriesReader< ImageType >        ReaderType;//读序列图片
//    ReaderType::Pointer itkReader = ReaderType::New();
//    typedef itk::GDCMImageIO       ImageIOType;//读DICOM图片
//    ImageIOType::Pointer dicomIO = ImageIOType::New();
//    itkReader->SetImageIO( dicomIO );//数据读入内存
//    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
//    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
//    nameGenerator->SetUseSeriesDetails( true );
//    nameGenerator->AddSeriesRestriction("0008|0021" );
//    nameGenerator->SetDirectory(filename);//设置文件目录
//
//    typedef std::vector< std::string >    SeriesIdContainer;
//    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
//    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();//迭代器
//    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
//    std::string seriesIdentifier;
//    seriesIdentifier = seriesUID.begin()->c_str();//通过迭代器读取所有单张切片
//    typedef std::vector< std::string >   FileNamesContainer;
//    FileNamesContainer fileNames;
//    fileNames = nameGenerator->GetFileNames( seriesIdentifier );
//    itkReader->SetFileNames( fileNames );
//    itkReader->Update();
//    dicomIO->GetMetaDataDictionary();//获取DIOCM头文件中信息
//
//    image = itkReader->GetOutput();
//    return EXIT_SUCCESS;
//}

int readITKImage( const string filename, ImageType::Pointer &image ){



    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer itkReader = ReaderType::New();
    itkReader->SetFileName( filename );
    itkReader->Update();
    image = itkReader->GetOutput();
    return EXIT_SUCCESS;
}


int main(int argc,char *argv[])
{
    string imageFileNameSource=string(argv[1]);
    string imageFileNameTarget=string(argv[2]);
    string outpath="/data/DDSM/normalization/benign/";
    vector<string> filename;
    boost::split( filename, imageFileNameSource, boost::is_any_of( "\\/" ), boost::token_compress_on );//如果此处使用boost::token_compress_off
    string name;
    name = filename[5];
    string outname = outpath+name;
    std::cout<<"filename "<<outname<<std::endl;

    ImageType::Pointer inputImageShortSource;
    readITKImage( imageFileNameSource,inputImageShortSource);

    ImageType::Pointer inputImageShortTarget;
    readITKImage( imageFileNameTarget,inputImageShortTarget);


    ImageType::SizeType size = inputImageShortSource->GetLargestPossibleRegion().GetSize();

    std::cout<<"image size is "<<size<<std::endl;

    matchingTwoImages(outname,inputImageShortSource, inputImageShortTarget);

}

