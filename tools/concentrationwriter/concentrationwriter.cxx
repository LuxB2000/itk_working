/*
 * constructors
 */
template<class TInputImage, const unsigned int Dimension>
ConcentrationWriter<TInputImage,Dimension>::ConcentrationWriter( ) {
	
}

template<class TInputImage, const unsigned int Dimension>
ConcentrationWriter<TInputImage,Dimension>::ConcentrationWriter(const char * path ) {
	m_path = std::string( path );
	ConcentrationWriter();
}

/*
 * destructor
 */
template<class TInputImage, const unsigned int Dimension>
ConcentrationWriter<TInputImage,Dimension>::~ConcentrationWriter( ) {
}

/*
 * SetInput
 */
template<class TInputImage, const unsigned int Dimension>
void
ConcentrationWriter<TInputImage,Dimension>::SetInput(InputImagePointerType im ) {
	m_image = im;
}

/*
 * Update
 */
template<class TInputImage, const unsigned int Dimension>
void
ConcentrationWriter<TInputImage,Dimension>::Update( ) {
		typedef itk::CastImageFilter< InputImageType, OutputImageType >
			CastFilterType;
		typename CastFilterType::Pointer cast = CastFilterType::New();
		cast->SetInput( m_image );

		typename MetaWriterType::Pointer writer = MetaWriterType::New();
		writer->SetFileName( m_path );
		writer->SetInput( cast->GetOutput() );
		writer->Update();
}



