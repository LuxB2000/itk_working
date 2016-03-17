/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description: Read an input 4D Data and extract statistics.
 *								Statistics are saved in a file or display on the console.
 *								The regions to compute the stat. are binary maps specified
 *								by the user by given the path to reach the folder where 
 *								they are saved.
 *
 *        Version:  1.0
 *        Created:  21/10/14 18:02:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

// my file readers
#include "../dcmwriter/dcmwriter.h"
#include "../dcmreader/dcmreader.h"
#include "../volume4dreader/volume4dreader.h"
#include "../console_tools/color.h"

#include "../meshtools/meshtools.h"

// general
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageIteratorWithIndex.h"
#include "itkImageConstIteratorWithIndex.h"

// folder parser
#include <sys/types.h>
#include <dirent.h>
#include <vector>
#include <math.h>

// mesh
#include "itkMesh.h"
#include "itkMeshFileReader.h"
#include "itkMeshSpatialObject.h"
#include <iostream>


/*
 * main
 */
int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 3 )
    {
    error( "Missing Parameters" );
    std::cerr << "Usage:" << argv[0]
			<< " input4DImagePath MapsFolder [removeOutliers, default:1] [outputFile, defautl: console]"
			<<std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvInput4dData = 1, argvInputFolder = 2, argcOutputFile=4, argcOutliers = 3;

	std::cout << "Input folder to parse: " << argv[argvInputFolder] 
		<< std::endl;

	bool removing_outlayers = true;
	const int COEF = 2;
	if( argc == 5 ){
		removing_outlayers = (atoi(argv[argcOutliers]) == 1);
	}

	if( removing_outlayers ){
		std::string m = std::string( "Removing outliers (out of born mean+/-" ) + SSTR(COEF) + std::string("*std)");
		blueMessage( m );
	}else{
		blueMessage("Doesn't remove outliers");
	}

	std::ofstream fileStream;
	bool writingFile = false;
	std::string outputFileName;
	if( argc>argcOutputFile ){
		outputFileName = std::string( argv[argcOutputFile] );
		writingFile = true;
		// write the information at the end of 
		fileStream.open(outputFileName.c_str(), std::ofstream::app );
		if( !fileStream.is_open() ){
			std::cerr << "Impossible to write statistics in file: "
				<< outputFileName << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	if( writingFile ){
		std::cout << "Print results in " << outputFileName << std::endl;
	}else{
		std::cout << "Results displayed in console." << std::endl;
	}
	std::cout << std::endl;


	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
	typedef  float          ArrayValType;
  typedef  float          PixelType;
  typedef  float          MeshPixelType;

	// image types
	typedef itk::VariableLengthVector<ArrayValType> ArrayType;
	typedef itk::VectorImage<ArrayValType,Dimension> ImageArrayType;
	ImageArrayType::Pointer volume4D;
	ImageArrayType::RegionType region;

	typedef itk::ImageRegionIterator< ImageArrayType > IteratorArrayType;
  typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::ImageRegionIterator< ImageArrayType > IteratorType;
	typedef itk::ImageIteratorWithIndex< ImageArrayType >
																													IteratorIndexType;
	ImageType::Pointer image;
	
	typedef DcmReader< ImageType, Dimension > DcmReaderType;
	DcmReaderType reader;

	// some data
	std::vector< PixelType > roiValues = std::vector<PixelType>();

	// compute statistics along each ROIs
	typedef std::vector<double> VectorDType;
	typedef std::vector<VectorDType> VectorVDType;

	//
	// Read the 4D data
	//
	typedef Volume4dReader< ImageArrayType, Dimension > VolumeReaderType;
	VolumeReaderType reader4D = VolumeReaderType();
	reader4D.SetFileName( argv[argvInput4dData] );
	try{
		reader4D.Update();
		volume4D = reader4D.GetOutput();
	}catch(itk::ExceptionObject &err ){
		std::cerr << "Error while reading the input 4D volume with path: "
			<< argv[argvInput4dData] << std::endl;
		std::cerr << err << std::endl;
	}
	int nbrOfComponent = volume4D->GetNumberOfComponentsPerPixel();

	//
	// parse folder and find all "mha" volume
	//

	struct dirent *de=NULL;
  DIR *d=NULL;
	d=opendir(argv[argvInputFolder]);
  if(d == NULL)
  {
		std::cerr << "Couldn't open input directory " <<
			argv[argvInputFolder] << std::endl;
    return(EXIT_FAILURE);
  }
	// Loop while not NULL
	std::string fileName, ext = std::string("mha");
	std::vector<std::string> fileToRead;
  while(de = readdir(d)){
		fileName= de->d_name;
		if( fileName.substr(fileName.find_last_of(".") +1 ) == ext ){
			fileToRead.push_back( fileName );
		}
	}

	const unsigned int nmbrOfFile = fileToRead.size();

	std::sort( fileToRead.rbegin(),
			fileToRead.rend(),
			std::greater<std::string>());

	std::cout << "Input binary masks: ";
	for( int i=0; i<fileToRead.size() ; i++){
		if( i!=0 ) std::cout << ", ";
		std::cout << fileToRead.at(i);
	}
	std::cout << std::endl;


	for( unsigned int f=0; f<nmbrOfFile; f ++ ){
		//
		// Read each binary mask
		//
		std::cout << std::endl;
		std::string path = std::string( std::string(argv[argvInputFolder]) +
				"/" + fileToRead.at(f) );
		reader = DcmReaderType();
		reader.SetFileName( path.c_str() );
		try{
			reader.Update();
			std::cout << "Reading file: " << path << std::endl;
			image = reader.GetOutput();
			image->DisconnectPipeline();
		}
		catch( itk::ExceptionObject & err ) {
			exit(EXIT_FAILURE);
		}

		// reset the roiValues
		roiValues = std::vector<PixelType>();
		
		//
		// Find how many regions are in the atlas
		//
		ImageArrayType::IndexType ind;
		ImageArrayType::PixelType vectConc;
		ImageType::PixelType roiVal;
		double res[2];
		typedef itk::ImageRegionConstIteratorWithIndex< ImageType > 
																	ConstIteratorType;
		// iterator parsing the BM
		ConstIteratorType it( image, image->GetLargestPossibleRegion() );
		for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
			if( it.Get() != 0 ){
				// check if the vector already contains the value
				if((std::find(roiValues.begin(), roiValues.end(),it.Get()) == 
						roiValues.end()) ) {
					roiValues.push_back( it.Get() );
				}
			}
		}

		std::cout << "We find the following ROIs: ";
		for( std::vector<PixelType>::const_iterator i = roiValues.begin();
				i != roiValues.end(); ++i){
			std::cout << *i << ' ';
		}
		std::cout << std::endl;


		VectorVDType mean = VectorVDType();
		VectorVDType N = VectorVDType();
		VectorVDType std  = VectorVDType();
		VectorVDType max_  = VectorVDType();
		VectorVDType min_  = VectorVDType();

		for( unsigned int i=0; i<roiValues.size(); i++ ){
			VectorDType vm = VectorDType(nbrOfComponent,0);
			mean.push_back( vm );
			VectorDType vs = VectorDType(nbrOfComponent,0);
			std.push_back( vs );
			VectorDType vn = VectorDType(nbrOfComponent,0);
			N.push_back( vn );
			VectorDType vmax = VectorDType(nbrOfComponent,0);
			max_.push_back( vmax );
			VectorDType vmin = VectorDType(nbrOfComponent,1000);
			min_.push_back( vmin );
		}

		//
		// measure the mean
		//
		int c=0;
		for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
			// at site in the region
			if( it.Get() != 0 ){
				for( c=0; c<roiValues.size(); c++ ){
					// find the correct vector place
					if(roiValues.at(c)==it.Get()){
						//VectorDType v = VectorDType(mean.at(c));
						//VectorDType n = VectorDType(N.at(c));
						vectConc = volume4D->GetPixel( it.GetIndex() );
						for(unsigned int i=0; i<nbrOfComponent; i++){
							if( vectConc[i] > 0 ){ // 1 ){
								mean.at(c).at(i) += vectConc[i];
								N.at(c).at(i) += 1;

								if( vectConc[i] < min_.at(c).at(i) ){
									min_.at(c).at(i) = vectConc[i];
								}
								if( vectConc[i] > max_.at(c).at(i) ){
									max_.at(c).at(i) = vectConc[i];
								}
								

							}
						}
						//mean.at(c) = VectorDType( v );
						//N.at(c) = VectorDType( n );
					}
				}
			}
		}

		for(c=0; c<roiValues.size(); c++){
			//VectorDType v = VectorDType( mean.at(c) );
			//VectorDType n = VectorDType( N.at(c) ) ;
			for(unsigned int i=0; i<nbrOfComponent; i++){
				//v.at(i) = v.at(i)/n.at(i);
					mean.at(c).at(i) = mean.at(c).at(i) / N.at(c).at(i);
			}
			//mean.at(c) = VectorDType( v );
		}

		//
		// compute std
		//
		for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
			if( it.Get() != 0 ){
				for( c=0; c<roiValues.size(); c++ ){
					if(roiValues.at(c)==it.Get()){
						//VectorDType v = std.at(c);
						//VectorDType m = mean.at(c);
						vectConc = volume4D->GetPixel( it.GetIndex() );
						for(unsigned int i=0; i<nbrOfComponent; i++){
							if( vectConc[i] > 0 ){ // 1 ){
								std.at(c).at(i) += pow(vectConc[i]-mean.at(c).at(i),2.0);
							}
						}
						//std.at(c) = v;
					}
				}
			}
		}
		for(c=0; c<roiValues.size(); c++){
			for(unsigned int i=0; i<nbrOfComponent; i++){
				std.at(c).at(i) = sqrt( std.at(c).at(i)/N.at(c).at(i) );
			}
		}
		

		VectorVDType mean_to_write = VectorVDType();
		VectorVDType std_to_write = VectorVDType();
		VectorVDType max_to_write  = VectorVDType();
		N = VectorVDType();
		VectorVDType min_to_write  = VectorVDType();
		if( removing_outlayers ){
			for( unsigned int i=0; i<roiValues.size(); i++ ){
				VectorDType vm = VectorDType(nbrOfComponent,0);
				mean_to_write.push_back( vm );
				VectorDType vs = VectorDType(nbrOfComponent,0);
				std_to_write.push_back( vs );
				VectorDType vn = VectorDType(nbrOfComponent,0);
				N.push_back( vn );
				VectorDType vmax = VectorDType(nbrOfComponent,0);
				max_to_write.push_back( vmax );
				VectorDType vmin = VectorDType(nbrOfComponent,1000);
				min_to_write.push_back( vmin );
			}

			// recompute all the stats: ignoring all the data that are outside mean+/2*std
			//
			// measure the mean
			//
			int c=0;
			for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
				// at site in the region
				if( it.Get() != 0 ){
					for( c=0; c<roiValues.size(); c++ ){
						// find the correct vector place
						if(roiValues.at(c)==it.Get()){
							vectConc = volume4D->GetPixel( it.GetIndex() );
							for(unsigned int i=0; i<nbrOfComponent; i++){
								if( (vectConc[i] < (mean.at(c).at(i)+COEF*std.at(c).at(i))) && (vectConc[i]>(mean.at(c).at(i))-COEF*std.at(c).at(i))  ){
									mean_to_write.at(c).at(i) += vectConc[i];
									N.at(c).at(i) += 1;

									if( vectConc[i] < min_to_write.at(c).at(i) ){
										min_to_write.at(c).at(i) = vectConc[i];
									}
									if( vectConc[i] > max_to_write.at(c).at(i) ){
										max_to_write.at(c).at(i) = vectConc[i];
									}
									

								}//end if in range mean+/-2*std
							}
						}
					}
				}
			}

			for(c=0; c<roiValues.size(); c++){
				for(unsigned int i=0; i<nbrOfComponent; i++){
						mean_to_write.at(c).at(i) = mean_to_write.at(c).at(i) / N.at(c).at(i);
				}
			}

			//
			// compute std
			//
			for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
				if( it.Get() != 0 ){
					for( c=0; c<roiValues.size(); c++ ){
						if(roiValues.at(c)==it.Get()){
							vectConc = volume4D->GetPixel( it.GetIndex() );
							for(unsigned int i=0; i<nbrOfComponent; i++){
								if( (vectConc[i] < (mean.at(c).at(i)+COEF*std.at(c).at(i))) && (vectConc[i]>(mean.at(c).at(i))-COEF*std.at(c).at(i))  ){
									std_to_write.at(c).at(i) += pow(vectConc[i]-mean_to_write.at(c).at(i),2.0);
								}
							}
						}
					}
				}
			}
			for(c=0; c<roiValues.size(); c++){
				for(unsigned int i=0; i<nbrOfComponent; i++){
					std_to_write.at(c).at(i) = sqrt( std_to_write.at(c).at(i)/N.at(c).at(i) );
				}
			}
		}else{ // we don't remove the outliers
			mean_to_write = mean;
			std_to_write = std;
			max_to_write = max_;
			min_to_write = min_;
		}
	
		if( !writingFile ){
			std::cout << "Statistics: [MeanTime1],[StdTime1],[min],[max];[MeanTime2],\n";
			for(c=0; c<roiValues.size(); c++){
				std::cout << " - " << roiValues.at(c) << ":";
				for( unsigned int i=0; i<nbrOfComponent; i++ ){
					std::cout << mean_to_write.at(c).at(i) << ", "
						" " << std_to_write.at(c).at(i) << ","
						" " << min_to_write.at(c).at(i)<<", " <<
						" " << max_to_write.at(c).at(i)<<", "
						<< ";";
				}
				std::cout << "\n";
			}
		}else{
			for(c=0; c<roiValues.size(); c++){
				fileStream << roiValues.at(c) << "\n";
				fileStream << "mean, std, min, max \n";
				for( unsigned int i=0; i<nbrOfComponent; i++ ){
					fileStream << mean_to_write.at(c).at(i) << ","
						" " << std_to_write.at(c).at(i) <<","
						" " << min_to_write.at(c).at(i) << ", " << max_to_write.at(c).at(i)
						<< "\n";
				}
				fileStream << "\n";
			}
			// write a separator for each file
			fileStream << " =======END OF FILE===========\n" ;
		}

	} // end for each mesh

	if( writingFile ){
		fileStream << "--------------------------------------------------------";
		fileStream << "--------------------------------------------------------\n";
		fileStream.close();
	}

	finalMessage( "All stats done.");

	std::cout << std::endl;
}

