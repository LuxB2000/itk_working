/*
 * constructors
 */
template<class TInputImage, const unsigned int Dimension>
Volume4dWriter<TInputImage,Dimension>::Volume4dWriter( ) {
	
}

template<class TInputImage, const unsigned int Dimension>
Volume4dWriter<TInputImage,Dimension>::Volume4dWriter(const char * path ) {
	m_path = std::string( path );
	Volume4dWriter();
}

/*
 * destructor
 */
template<class TInputImage, const unsigned int Dimension>
Volume4dWriter<TInputImage,Dimension>::~Volume4dWriter( ) {
}

/*
 * SetInput
 */
template<class TInputImage, const unsigned int Dimension>
void
Volume4dWriter<TInputImage,Dimension>::SetInput(InputImagePointerType im ) {
	m_image = im;
}

/*
 * Update
 */
template<class TInputImage, const unsigned int Dimension>
void
Volume4dWriter<TInputImage,Dimension>::Update( ) {
		typedef itk::CastImageFilter< InputImageType, OutputImageType >
			CastFilterType;
		typename CastFilterType::Pointer cast = CastFilterType::New();
		cast->SetInput( m_image );

		typename MetaWriterType::Pointer writer = MetaWriterType::New();
		writer->SetFileName( m_path );
		writer->SetInput( cast->GetOutput() );
		writer->Update();
}



