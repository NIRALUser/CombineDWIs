#include "CombineDWIsCLP.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageAlgorithm.h>
#include <itkMetaDataObject.h>
#include <cmath>

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1


// Separate the vector image into a vector of images
template <class PixelType>
int SeparateImages( const typename itk::VectorImage<PixelType, 3>
                    ::Pointer & imagePile,
                    std::vector<typename itk::Image<PixelType, 3>::Pointer> & vectorImage
                    )
{
    typedef itk::Image<PixelType, 3> ImageType;
    typedef itk::VectorImage<PixelType, 3> VectorImageType;
    typename itk::VectorImage<PixelType, 3>::SizeType size;
    typename itk::VectorImage<PixelType, 3>::DirectionType direction;
    typename itk::VectorImage<PixelType, 3>::PointType origin;
    typename itk::VectorImage<PixelType, 3>::SpacingType spacing;
    size = imagePile->GetLargestPossibleRegion().GetSize();
    direction = imagePile->GetDirection();
    origin = imagePile->GetOrigin();
    spacing = imagePile->GetSpacing();
    typename itk::ImageRegionIterator<VectorImageType> in( imagePile,
                                                           imagePile->GetLargestPossibleRegion() );
    typedef typename itk::ImageRegionIterator<ImageType> IteratorImageType;
    std::vector<IteratorImageType> out;
    for( unsigned int i = 0; i < imagePile->GetVectorLength(); i++ )
    {
        typename ImageType::Pointer imageTemp = ImageType::New();
        imageTemp->SetRegions( size );
        imageTemp->SetOrigin( origin );
        imageTemp->SetDirection( direction );
        imageTemp->SetSpacing( spacing );
        imageTemp->Allocate();
        vectorImage.push_back( imageTemp );
        IteratorImageType outtemp( imageTemp, imageTemp->GetLargestPossibleRegion() );
        outtemp.GoToBegin();
        out.push_back( outtemp );
    }
    for( in.GoToBegin(); !in.IsAtEnd(); ++in )
    {
        itk::VariableLengthVector<PixelType> value = in.Get();
        for( unsigned int i = 0; i < imagePile->GetVectorLength(); i++ )
        {
            out[i].Set( value[i] );
            ++out[i];
        }
    }
    return EXIT_SUCCESS;
}


template<class PixelType>
void
ReadInputImage( std::string fileName ,
                std::vector< typename itk::Image< PixelType , 3 >::Pointer > &vectorOfImage ,
                itk::MetaDataDictionary &dict
               )
{
    // open image file
    typedef itk::VectorImage< PixelType , 3 > VectorImageType ;
    typename itk::ImageFileReader< VectorImageType >::Pointer reader ;
    reader = itk::ImageFileReader< VectorImageType >::New() ;
    reader->SetFileName( fileName.c_str() );
    reader->Update() ;
    // Save metadata dictionary
    dict = reader->GetOutput()->GetMetaDataDictionary();
    //Separate the vector image into a vector of images
    SeparateImages<PixelType>( reader->GetOutput(), vectorOfImage );
}

// Write back the vector of images into a image vector
template <class PixelType>
int AddImage( typename itk::VectorImage<PixelType, 3>
              ::Pointer & imagePile,
              const std::vector<typename itk::Image<PixelType, 3>::Pointer> & vectorImage
              )
{
    typedef itk::Image<PixelType, 3> ImageType;
    imagePile->SetRegions( vectorImage.at( 0 )->GetLargestPossibleRegion().GetSize() );
    imagePile->SetOrigin( vectorImage.at( 0 )->GetOrigin() );
    imagePile->SetDirection( vectorImage.at( 0 )->GetDirection() );
    imagePile->SetSpacing( vectorImage.at( 0 )->GetSpacing() );
    imagePile->SetVectorLength( vectorImage.size() );
    imagePile->Allocate();
    typename itk::ImageRegionIterator<itk::VectorImage<PixelType, 3> > out( imagePile,
                                                                            imagePile->GetLargestPossibleRegion()
                                                                            );
    typedef typename itk::ImageRegionIterator<ImageType> IteratorImageType;
    std::vector<IteratorImageType> in;
    for( unsigned int i = 0; i < imagePile->GetVectorLength(); i++ )
    {
        IteratorImageType intemp( vectorImage.at( i ), vectorImage.at( i )->GetLargestPossibleRegion() );
        intemp.GoToBegin();
        in.push_back( intemp );
    }
    itk::VariableLengthVector<PixelType> value;
    value.SetSize( vectorImage.size() );
    for( out.GoToBegin(); !out.IsAtEnd(); ++out )
    {
        for( unsigned int i = 0; i < imagePile->GetVectorLength(); i++ )
        {
            value.SetElement( i, in.at( i ).Get() );
            ++in[i];
        }
        out.Set( value );
    }
    return EXIT_SUCCESS;
}

int GetBValueFromMetaDataDictionary( const itk::MetaDataDictionary &dict , double &bval )
{
    std::string bvalueString ;
    itk::ExposeMetaData<std::string>( dict , "DWMRI_b-value" , bvalueString ) ;
    try
    {
      bval = std::stod( bvalueString ) ;
    }
    catch(...)
    {
        std::cerr << "Cannot convert b-value to double" << std::endl ;
        std::cerr << "b-value: " << bvalueString << std::endl ;
        return 1 ;
    }
    return 0 ;
}

int FindLastGradient( const itk::MetaDataDictionary &dict )
{
  typedef itk::MetaDataObject<std::string> MetaDataStringType;
  itk::MetaDataDictionary::ConstIterator itr = dict.Begin();
  itk::MetaDataDictionary::ConstIterator end = dict.End();
  int maxVal = -1 ;
  while( itr != end )
  {
      itk::MetaDataObjectBase::Pointer entry = itr->second;
      MetaDataStringType::Pointer entryvalue
              = dynamic_cast<MetaDataStringType *>( entry.GetPointer() );
      if( entryvalue )
      {
          // get the gradient directions
          int pos = itr->first.find( "DWMRI_gradient" );
          if( pos != -1 )
          {
              std::string gradientNumberString =  itr->first.substr( 15 ) ; //15 characters in "DWMRI_gradient_"
              int gradientNumber = std::stoi( gradientNumberString ) ;
              if( gradientNumber > maxVal )
              {
                  maxVal = gradientNumber ;
              }
          }
      }
      ++itr;
  }
  return maxVal ;
}

std::vector< itk::Vector< double, 3> > ExtractGradientsFromMetaDataDictionary( const itk::MetaDataDictionary &dict )
{
    std::vector< itk::Vector< double, 3> > gradients ;
    typedef itk::MetaDataObject<std::string> MetaDataStringType;
    itk::MetaDataDictionary::ConstIterator itr = dict.Begin();
    itk::MetaDataDictionary::ConstIterator end = dict.End();
    while( itr != end )
    {
        itk::MetaDataObjectBase::Pointer entry = itr->second;
        MetaDataStringType::Pointer entryvalue
                = dynamic_cast<MetaDataStringType *>( entry.GetPointer() );
        if( entryvalue )
        {
            // get the gradient directions
            int pos = itr->first.find( "DWMRI_gradient" );
            if( pos != -1 )
            {
                std::string tagvalue = entryvalue->GetMetaDataObjectValue();
                itk::Vector<double, 3> vec;
                std::istringstream iss( tagvalue );
                iss >> vec[0] >> vec[1] >> vec[2]; // we copy the metavalue in an itk::vector
                if( iss.fail() )
                {
                    iss.str( tagvalue );
                    iss.clear();
                    std::string trash;
                    iss >> vec[0] >> trash >> vec[1] >> trash >> vec[2]; // in case the separator between the values is
                    // something else than spaces
                    if( iss.fail() ) // problem reading the gradient values
                    {
                        std::cerr << "Error reading a DWMRI gradient value - skipping" << itr->first << std::endl ;
                    }
                }
              gradients.push_back( vec ) ;
            }
        }
        ++itr;
    }
    return gradients ;
}

void
AppendDictionary( itk::MetaDataDictionary &dictAll , const itk::MetaDataDictionary &currentDict , bool verbose )
{
    double bvalueRef ;
    GetBValueFromMetaDataDictionary( dictAll , bvalueRef ) ;
    double currentbval ;
    GetBValueFromMetaDataDictionary( currentDict , currentbval ) ;
    double ratio = std::sqrt( currentbval / bvalueRef ) ;//multiply each component by sqrt of ratio!
    if( verbose )
    {
        std::cout << "Reference b-value: " << bvalueRef << std::endl ;
        std::cout << "Current image b-value: " << currentbval <<std::endl ;
        std::cout << "Ratio (sqrt):" << ratio << std::endl ;
    }
    int lastGradient = FindLastGradient( dictAll ) ;
    if( verbose )
    {
        std::cout << "Last gradient in appended dictionary: " << lastGradient << std::endl ;
    }
    std::vector< itk::Vector< double, 3> > currentGradientsVec ;
    currentGradientsVec = ExtractGradientsFromMetaDataDictionary( currentDict ) ;
    for( size_t i = 0 ; i < currentGradientsVec.size() ; i++ )
    {
        std::ostringstream oss;
        if( verbose )
        {
            std::cout << "Current image gradient " << i << std::endl ;
            std::cout << currentGradientsVec[ i ][ 0 ] << " " << currentGradientsVec[ i ][ 1 ] << " " <<  currentGradientsVec[ i ][ 2 ] << std::endl ;
        }
        // write the new gradient values (after transformation) in the metadatadictionary
        if( currentGradientsVec[i].GetNorm() > .00001 ) //we do not multiply by the ratio if it is a b0 image
        {
            for( int j = 0 ; j < 3 ; j++ )
            {
                currentGradientsVec[i][j] *= ratio ;
            }
        }
        if( verbose )
        {
            std::cout << "Current image gradient (with ratio applied): " << i << std::endl ;
            std::cout << currentGradientsVec[ i ][ 0 ] << " " << currentGradientsVec[ i ][ 1 ] << " " <<  currentGradientsVec[ i ][ 2 ] << std::endl ;
        }
        oss << currentGradientsVec[i][0] << " "
            << currentGradientsVec[i][1] << " "
            << currentGradientsVec[i][2] ;
        unsigned long long currentGradientIndex = lastGradient + 1 + i ;
        std::ostringstream sstagkey ;
        sstagkey << std::setw( 4 ) << std::setfill( '0' ) << currentGradientIndex;
        std::string tagkey = "DWMRI_gradient_" + sstagkey.str() ;
        if( verbose )
        {
            std::cout << "Adding following key to dictionary: " << tagkey << ":=" << oss.str() << std::endl ;
        }
        itk::EncapsulateMetaData<std::string>( dictAll , tagkey, oss.str() ) ;
    }
}


template<class PixelType>
int DoIt( int argc, char * argv[] )
{

    /////////////////   Get arguments using GenerateCLP parser to get all the arguments  ///////////////////////
    PARSE_ARGS ;
    typedef itk::Image<PixelType, 3> ImageType;
    typedef itk::VectorImage<PixelType, 3> VectorImageType;
    std::vector<typename ImageType::Pointer> vectorOfImages;
    itk::MetaDataDictionary dictAll;
    if( verbose )
    {
        std::cout << "Number of input volumes: " << inputVolumes.size() << std::endl ;
        std::cout << "Reading reference volume" << std::endl ;
        std::cout << "Reference Volume: " << inputVolumes[ 0] << std::endl ;
    }
    ReadInputImage<PixelType>( inputVolumes[0] , vectorOfImages , dictAll ) ;
    if( verbose )
    {
        std::cout << "Number of Images in reference DWI: " << vectorOfImages.size() << std::endl ;
    }
    for( unsigned int i = 1 ; i < inputVolumes.size() ; i++ )
    {
        if( verbose )
        {
            std::cout << "Processing volume: " << inputVolumes[ i ] << std::endl ;
        }
       std::vector<typename ImageType::Pointer> currentVectorOfImages;
       itk::MetaDataDictionary currentDict;
       ReadInputImage<PixelType>( inputVolumes[ i ] , currentVectorOfImages , currentDict ) ;
       vectorOfImages.insert( vectorOfImages.end(), currentVectorOfImages.begin(), currentVectorOfImages.end()) ;
       AppendDictionary( dictAll , currentDict , verbose ) ;
    }
    typename itk::VectorImage< PixelType, 3 >::Pointer outputImage ;
    outputImage = itk::VectorImage<PixelType, 3>::New() ;
    if( verbose )
    {
      std::cout << "Converting std::vector of itk::Images, to itk::VectorImage" << std::endl ;
    }
    AddImage<PixelType>( outputImage, vectorOfImages ) ;
    outputImage->SetMetaDataDictionary( dictAll );
    vectorOfImages.clear() ;
    if( verbose )
    {
        std::cout << "Output Volume name: " << outputVolume << std::endl ;
        std::cout << "Writing output volume" << std::endl ;
    }
    typedef itk::ImageFileWriter<VectorImageType> WriterType ;
    try
    {
        typename WriterType::Pointer writer = WriterType::New() ;
        writer->SetInput( outputImage ) ;
        writer->SetFileName( outputVolume.c_str() ) ;
        writer->UseCompressionOn() ;
        writer->Update() ;
    }
    catch( itk::ExceptionObject exception )
    {
        std::cerr << exception << std::endl ;
        return EXIT_FAILURE ;
    }
    return EXIT_SUCCESS ;
}

void GetImageType( std::string fileName,
                   itk::ImageIOBase::IOPixelType &pixelType ,
                   itk::ImageIOBase::IOComponentType &componentType
                 )
{
  typedef itk::Image< unsigned char , 3 > ImageType ;
  itk::ImageFileReader< ImageType >::Pointer imageReader 
              = itk::ImageFileReader< ImageType >::New() ;
  imageReader->SetFileName( fileName.c_str() ) ;
  imageReader->UpdateOutputInformation() ;
  pixelType = imageReader->GetImageIO()->GetPixelType() ;
  componentType = imageReader->GetImageIO()->GetComponentType() ;
}



int main( int argc, char * argv[] )
{
  /////////////////   Get arguments using GenerateCLP
  //// parser to get the input volume filename  ///////////////////////
  PARSE_ARGS;
  if( inputVolumes.empty() )
  {
    std::cerr << "At least one input image need to be specified" << std::endl ;
    return EXIT_FAILURE ;
  }
  if( outputVolume.empty() )
  {
    std::cerr << "An output file name needs to be specified" << std::endl ;
    return EXIT_FAILURE ;
  }
  /////////////////   Read input volume data type and instantiate
  //// the corresponding templated filter  function    ////////////////////
  itk::ImageIOBase::IOPixelType pixelType ;
  itk::ImageIOBase::IOComponentType componentType ;
  try
  {
    GetImageType ( inputVolumes[0] , pixelType , componentType ) ;
    // This filter handles all image component types
    switch( componentType )
    {
      case itk::ImageIOBase::UCHAR:
        return DoIt< unsigned char >( argc , argv ) ;
      case itk::ImageIOBase::CHAR:
        return DoIt< char >( argc , argv ) ;
      case itk::ImageIOBase::USHORT:
        return DoIt< unsigned short >( argc , argv ) ;
      case itk::ImageIOBase::SHORT:
        return DoIt< short >( argc , argv ) ;
      case itk::ImageIOBase::UINT:
        return DoIt< unsigned int >( argc , argv ) ;
      case itk::ImageIOBase::INT:
        return DoIt< int >( argc , argv ) ;
      case itk::ImageIOBase::ULONG:
        return DoIt< unsigned long >( argc , argv ) ;
      case itk::ImageIOBase::LONG:
        return DoIt< long >( argc , argv ) ;
      case itk::ImageIOBase::FLOAT:
        return DoIt< float >( argc , argv ) ;
      case itk::ImageIOBase::DOUBLE:
        return DoIt< double >( argc , argv ) ;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
        std::cout << "unknown component type" << std::endl ;
        return EXIT_FAILURE ;
        break ;
    }
  }
  catch( itk::ExceptionObject &excep )
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl ;
    std::cerr << excep << std::endl ;
    return EXIT_FAILURE ;
  }
  //should never arrive here
  return EXIT_FAILURE ;
}
