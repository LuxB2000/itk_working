#include "dcmreader.h"

/*
 * constructors
 */
template<class TOutputImage, const unsigned int Dimension>
DcmReader<TOutputImage,Dimension>::DcmReader( ) {
	m_param_path = std::string();
}

template<class TOutputImage, const unsigned int Dimension>
DcmReader<TOutputImage,Dimension>::DcmReader( const char* p ) {
	m_path = std::string(p);
	DcmReader();
}

template<class TOutputImage, const unsigned int Dimension>
void
DcmReader<TOutputImage,Dimension>::m_parseName() {

	m_type=0;
	// type=0 no valid path
	// type=1 if MHA image (extension .mha or .mhd)
	// type=2 if 3D tiff image (extension .tif)
	// type=99 if path designates a folder

	std::vector<std::string> extensions = std::vector<std::string>();
	extensions.push_back( std::string(".mha") );
	extensions.push_back( std::string(".mhd") );
	extensions.push_back( std::string(".tif") );


	if( !m_path.empty() ){
		struct stat status;
		if( stat( m_path.c_str(), &status ) == -1 ){
			perror("stat");
			itk::ExceptionObject excp;
      excp.SetDescription("The path " + m_path + 
					" is not valide.");
      throw excp;
		}
		if( status.st_mode & S_IFDIR )
		{ // FOLDER
			m_type=99;
		}else if( status.st_mode & S_IFREG )
		{ // FILE
			//check extension
			std::string ext = m_path.substr(m_path.size()-4,m_path.size());
			//std::cout<<"TEST:-"<<ext<<"-"<<std::endl;
			if( (ext.compare( extensions[0] ) == 0) || 
				(ext.compare(extensions[1]) == 0) )
				{ // META image
					m_type=1;
				}
			if( ext.compare( extensions[2] ) == 0 )
			{ // tif image
				m_type=2;
			}
		}
	}
	std::cout<<"The path reaches";
	switch(m_type){
	case 1:
		std::cout<<" a META file";
		break;
	case 2:
		std::cout<<" a TIF file";
		break;
	case 99:
		std::cout<<" a folder";
		break;
	default:
		std::cout<<" NOTHING VALID";
		break;
	}
	std::cout<<std::endl;

	if( m_type == 0 )
	{
		itk::ExceptionObject excp;
		excp.SetDescription("The path " + m_path + 
				" does not reach any valid file");
		throw excp;
	}

	m_output = 0; 
}

/*
 * destructor
 */
template<class TOutputImage, const unsigned int Dimension>
DcmReader<TOutputImage,Dimension>::~DcmReader(){
}

/*
 * Update
 */
template<class TOutputImage, const unsigned int Dimension>
void
DcmReader<TOutputImage,Dimension>::Update(){
	
	m_parseName();

	if( m_type == 1 ) // META file
	{
		MetaReaderType::Pointer     reader = MetaReaderType::New();
		reader->SetFileName( m_path );
		// cast the image
		typedef itk::CastImageFilter< 
                        MetaImageType,
                        OutputImageType >     CastFilterType;
		typename CastFilterType::Pointer caster = CastFilterType::New();
		caster->SetInput( reader->GetOutput() );

		//try{
			caster->Update();
			m_output = caster->GetOutput();
		//}
		//catch( itk::ExceptionObject & err ) 
		//{
		//	m_output = 0;
		//	std::cerr << "Error while reading and casting the input Meta file!";
		//	std::cerr << std::endl;
		//	std::cerr << err << std::endl; 
		//}
	}
	else if(m_type == 2 ) // 3D TIF volumes
	{
		TIFFReaderType::Pointer reader = TIFFReaderType::New();
		reader->SetFileName( m_path );
		typedef itk::CastImageFilter< TIFFImageType,
                                  OutputImageType> CastFilterType;
		typename CastFilterType::Pointer caster = CastFilterType::New();
		caster->SetInput( reader->GetOutput() );
		//try{
			caster->Update();
			m_output = caster->GetOutput();
		//}
		//catch( itk::ExceptionObject & err ){
		//	m_output = 0;
		//	throw excp;
		//}

		//
		// extract parameters from TXT files
		//
		// FIXME path is extracted as follow: all the path BEFORE /3D-data/ but the
		// number '3' could be used a fileID .. so we use the last 'D'. It should be the
		// complete string '3D-data' ...
		std::string path = m_path.substr( 0,m_path.find_last_of(std::string("D"))-1 );
		std::string tifName = m_path.substr( m_path.find_last_of("/\\") +1 );
		std::string tifID = tifName.substr( 0,tifName.find_first_of("_") ); // WARNING: Make a change due to new data, the file is ID_-_PostPd_T1.tif
		// find the text file
		//std::cout << "test: " << path << " " << tifName << " " << tifID << std::endl;
		path += std::string("Parameters/") + tifID + "-parameters.txt";
		//std::cout<<"TEST:"<<path<<std::endl;
		struct stat status;
		typename OutputImageType::SpacingType spacing;
		typename OutputImageType::PointType origin;

		if( stat( path.c_str(), &status ) == -1 ){
			std::cerr << "We CAN NOT find the parameters file for the TIF image!\n";
			std::cerr << "Parameter file with path: " << path << "\n";
			std::cerr<<"WARNING: we use default parameters (spacing, origin, ..)";
			std::cerr<<std::endl;

			spacing[0] = 1;
			spacing[1] = 1;
			spacing[2] = 1;

			origin[0] = 0;
			origin[1] = 0;
			origin[2] = 0;
		}else{
			// read the file
			ReadTiffImageParameters parametersReader =
				ReadTiffImageParameters( path );
			parametersReader.Read();
			TiffParameters* p = parametersReader.GetParamters();
			m_param_path = path;
			
			spacing[0] = p->spacing[0];
			spacing[1] = p->spacing[1];
			spacing[2] = p->spacing[2];

			origin[0] = p->origin[0];
			origin[1] = p->origin[1];
			origin[2] = p->origin[2];
		}

		m_output->SetSpacing( spacing );
		m_output->SetOrigin( origin );
			
		//typename OutputImageType::DirectionType direction;
		//direction = m_output->GetDirection();
		//std::cout<<"DirectiOns: "<<direction<<std::endl;

	}

	else if(m_type == 99 ) // 2D DCM slices
	{
		DCMReaderType::Pointer reader = DCMReaderType::New();
		ImageIOType::Pointer        dicomIO = ImageIOType::New();
		NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

		reader->SetImageIO( dicomIO );
 
		nameGenerator->SetUseSeriesDetails( true );
		// restrict to the serie date
		nameGenerator->AddSeriesRestriction("0008|0021" );
		nameGenerator->SetDirectory( m_path );

		const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
		std::string seriesIdentifier = seriesUID.begin()->c_str();
		FileNamesContainer fileNames;
		fileNames = nameGenerator->GetFileNames( seriesIdentifier );
		reader->SetFileNames( fileNames );


		// cast the image
		typedef itk::CastImageFilter< 
                        DCMImageType,
                        OutputImageType >          CastFilterType;
		typename CastFilterType::Pointer caster = CastFilterType::New();
		caster->SetInput( reader->GetOutput() );

		try{
			caster->Update();
			m_output = caster->GetOutput();
		}
		catch( itk::ExceptionObject & err ) 
		{
			m_output = 0;
			throw( err );
		}
	}
}

/*
 * GetOutput
 */
template<class TOutputImage, const unsigned int Dimension>
typename DcmReader<TOutputImage,Dimension>::OutputPointerType
DcmReader<TOutputImage,Dimension>::GetOutput(){
	return m_output;
}
