/*
 * Constructor
 */
template<class TInputImage, const unsigned int Dimension>
DcmWriter<TInputImage,Dimension>::DcmWriter( ){
}

template<class TInputImage, const unsigned int Dimension>
DcmWriter<TInputImage,Dimension>::DcmWriter( const char* p ){
	m_path = std::string( p );
	DcmWriter();
}

/*
 * Destructor
 */
template<class TInputImage, const unsigned int Dimension>
DcmWriter<TInputImage,Dimension>::~DcmWriter( ){
}

/*
 * SetInputImage
 */
template<class TInputImage, const unsigned int Dimension>
void
DcmWriter<TInputImage,Dimension>::SetInput( 
		typename DcmWriter<TInputImage,Dimension>::InputImagePointerType image){

	m_image = image;
	m_mode = 1;
}


/*
 * Update
 */
template<class TInputImage, const unsigned int Dimension>
void
DcmWriter<TInputImage,Dimension>::Update( ){
	if( m_mode == 1)
	{
		typedef itk::CastImageFilter< InputImageType, OutputImageType >
			CastFilterType;
		typename CastFilterType::Pointer cast = CastFilterType::New();
		cast->SetInput( m_image );

		typename MetaWriterType::Pointer writer = MetaWriterType::New();
		writer->SetFileName( m_path );
		writer->SetInput( cast->GetOutput() );
		//try{
			writer->Update();
		//}catch( itk::ExceptionObject &err ){
		//	std::cerr << "Error while writing the image as: "
		//		<< m_path << std::endl;
		//}
	}else if( m_mode == 2 )
	{
		typedef itk::NumericSeriesFileNames OutputNamesGeneratorType;
		typename SeriesWriterType::DictionaryArrayType dico = GenerateDico();
		ImageIOType::Pointer gdcmIO = ImageIOType::New();
		gdcmIO->SetNumberOfDimensions( Dimension );
		for( unsigned int i=0; i<Dimension; i++ )
		{
			//gdcmIO->SetDimensions( i, m_image->GetRequestedRegion().GetSize()[i] );
			gdcmIO->SetOrigin( i, m_image->GetOrigin()[i] );
			gdcmIO->SetSpacing( i, m_image->GetSpacing()[i] );
		}

		// Make the output directory and generate the file names.
		itksys::SystemTools::MakeDirectory( m_path.c_str() );
		// Generate the file names
		OutputNamesGeneratorType::Pointer outputNames =
																				OutputNamesGeneratorType::New();
		std::string seriesFormat(m_path);
		seriesFormat = seriesFormat + "/" + "IM%d.dcm";
		outputNames->SetSeriesFormat (seriesFormat.c_str());
		outputNames->SetStartIndex(1);
		outputNames->SetEndIndex(1);//m_image->GetRequestedRegion().GetSize()[2]);

		typename SeriesWriterType::Pointer seriesWriter =
																								SeriesWriterType::New();
		seriesWriter->SetInput( m_image );
    seriesWriter->SetImageIO( gdcmIO ); 
		seriesWriter->SetFileNames( outputNames->GetFileNames() );
    seriesWriter->SetMetaDataDictionaryArray( &dico ); 

		seriesWriter->Update();
	}else
	{
		itk::ExceptionObject excp;
     excp.SetDescription("Mode is not well specified is should be 1 or 2");
      throw excp;
	}
}



/*
 * Generate Image IO
 * See http://www.itk.org/Wiki/ITK/Examples/DICOM/ResampleDICOM
 */
template<class TInputImage, const unsigned int Dimension>
typename DcmWriter<TInputImage,Dimension>::SeriesWriterType::DictionaryArrayType
DcmWriter<TInputImage,Dimension>::GenerateDico()
{
	typedef typename SeriesWriterType::DictionaryRawPointer DictionaryType;
	typename SeriesWriterType::DictionaryArrayType outputDico;
	typename InputImageType::SizeType size = 
																		m_image->GetRequestedRegion().GetSize();

	for( unsigned int i=0; i<size[2]; i++ )
	{
		// Create a new dictionary for this slice
		DictionaryType dict = new typename SeriesWriterType::DictionaryType;
		#if ITK_VERSION_MAJOR >= 4
		gdcm::UIDGenerator suid;
		std::string seriesUID = suid.Generate();
		gdcm::UIDGenerator fuid;
		std::string frameOfReferenceUID = fuid.Generate();
#else
		std::string seriesUID = gdcm::Util::CreateUniqueUID( gdcmIO->GetUIDPrefix());
		std::string frameOfReferenceUID = gdcm::Util::CreateUniqueUID( gdcmIO->GetUIDPrefix());
#endif
		std::string studyUID = "";
		std::string sopClassUID = "";
		itk::ExposeMetaData<std::string>(*dict, "0020|000d", studyUID);
		itk::ExposeMetaData<std::string>(*dict, "0008|0016", sopClassUID);
		// Set the UID's for the study, series, SOP  and frame of reference
    itk::EncapsulateMetaData<std::string>(*dict,"0020|000e", seriesUID);
    itk::EncapsulateMetaData<std::string>(*dict,"0020|0052", frameOfReferenceUID);
#if ITK_VERSION_MAJOR >= 4
    gdcm::UIDGenerator sopuid;
    std::string sopInstanceUID = sopuid.Generate();
#else
    std::string sopInstanceUID = gdcm::Util::CreateUniqueUID( gdcmIO->GetUIDPrefix());
#endif
    itk::EncapsulateMetaData<std::string>(*dict,"0008|0018", sopInstanceUID);
    itk::EncapsulateMetaData<std::string>(*dict,"0002|0003", sopInstanceUID);


		// Change fields that are slice specific
    itksys_ios::ostringstream value;
    value.str("");
    value << i + 1;
 
    // Image Number
    itk::EncapsulateMetaData<std::string>(*dict,"0020|0013", value.str());

		// Series Description - Append new description to current series
    // description
		value.str("");
    value << "Manually generate by J. Plumat,"
					<< " contact: j.plumat[at]auckland.ac.nz";
		// This is an long string and there is a 64 character limit in the
    // standard
    unsigned lengthDesc = value.str().length();
 
    std::string seriesDesc( value.str(), 0,
                            lengthDesc > 64 ? 64
                            : lengthDesc);
    itk::EncapsulateMetaData<std::string>(*dict,"0008|103e", seriesDesc);

		// Series Number
    value.str("");
    value << 1001;
    itk::EncapsulateMetaData<std::string>(*dict,"0020|0011", value.str());
 
    // Derivation Description - How this image was derived
    value.str("");
    value << "Generate by ITK.";
 
    lengthDesc = value.str().length();
    std::string derivationDesc( value.str(), 0,
                                lengthDesc > 1024 ? 1024
                                : lengthDesc);
    itk::EncapsulateMetaData<std::string>(*dict,"0008|2111", derivationDesc);
 
    // Image Position Patient: This is calculated by computing the
    // physical coordinate of the first pixel in each slice.
    typename OutputImageType::PointType position;
    typename OutputImageType::IndexType index;
    index[0] = 0;
    index[1] = 0;
    index[2] = i;
    m_image->TransformIndexToPhysicalPoint(index, position);

 
    value.str("");
    value << position[0] << "\\" << position[1] << "\\" << position[2];
    itk::EncapsulateMetaData<std::string>(*dict,"0020|0032", value.str());
    // Slice Location: For now, we store the z component of the Image
    // Position Patient.
    value.str("");
    value << position[2];
    itk::EncapsulateMetaData<std::string>(*dict,"0020|1041", value.str());

		// Save the dictionary
    outputDico.push_back(dict);

	}
	return outputDico;
}


