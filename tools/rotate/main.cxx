/*
 * Rotate an image.
 * Center of rotation is in the middle of the volume
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 6 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputImagePath";
		std::cerr << " Angle in degress along axis X Y and Z" <<std::endl;
    return EXIT_FAILURE;
  }
	const unsigned int argcInputImage = 1, argcOutputImage = 2,
							 argcAngleX = 3, argcAngleY = 4, argcAngleZ = 5;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argcInputImage] );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();// dcmReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// get the angle specify by user
	const double deg2rad = std::atan(1.0)/45.0;
	const double Xdeg = atof(argv[argcAngleX]),
		Ydeg = atof(argv[argcAngleY]),
		Zdeg = atof(argv[argcAngleZ]);

	std::string msg = "";
	msg += "input angle degrees along axis X=" + SSTR(Xdeg)
		+ " Y=" + SSTR(Ydeg) + " Z= " + SSTR(Zdeg);
	blueMessage(msg);

	itk::Vector<double,3> x_axis, y_axis, z_axis;
	x_axis[0]=1;x_axis[1]=0;x_axis[2]=0;
	y_axis[0]=0;y_axis[1]=1;y_axis[2]=0;
	z_axis[0]=0;z_axis[1]=0;z_axis[2]=1;

	// set the transform parameters
	typedef itk::AffineTransform<double, Dimension> TransformType;
	const double Xrad = Xdeg * deg2rad;
	const double Yrad = Ydeg * deg2rad;
	const double Zrad = Zdeg * deg2rad;

	TransformType::Pointer transform = TransformType::New();

	// set the center of rotation
	typedef TransformType::CenterType CenterType;
	ImageType::IndexType ind;
	CenterType center;
	for(unsigned int i=0; i<2; i++){
		ind[i] =
			static_cast<unsigned int>(
				image->GetLargestPossibleRegion().GetSize()[i] / 2
				);
	}
	ind[2] = 0;
	image->TransformIndexToPhysicalPoint(ind,center);
	msg = "center of rotation: " + SSTR(center[0]) + " " 
		+ SSTR(center[1]) + " " + SSTR(center[2]) + "\n";
	blueMessage(msg);

	transform->SetCenter( center );

	transform->Rotate3D(x_axis,Xrad);
	transform->Rotate3D(y_axis,Yrad);
	transform->Rotate3D(z_axis,Zrad);

	typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
	ResampleImageFilterType::Pointer filter = ResampleImageFilterType::New();
	filter->SetInput(image);
	filter->SetSize(image->GetLargestPossibleRegion().GetSize());
	filter->SetOutputSpacing(image->GetSpacing());
	filter->SetOutputOrigin(image->GetOrigin());
	filter->SetOutputDirection(image->GetDirection());
	filter->SetTransform(transform);
	filter->UpdateLargestPossibleRegion();

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argcOutputImage] );
	writer.SetInput( filter->GetOutput() );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< argv[argcOutputImage] << std::endl;
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}
	
	msg = "Output volume wrote with ";
	msg += argv[argcOutputImage];
	finalMessage( msg );
	std::cout << std::endl;
}

