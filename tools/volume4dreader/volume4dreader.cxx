/*
 * constructors
 */
template<class TOutputImage, const unsigned int Dimension>
Volume4dReader<TOutputImage,Dimension>::Volume4dReader( ) {
	
}

template<class TOutputImage, const unsigned int Dimension>
Volume4dReader<TOutputImage,Dimension>::Volume4dReader(const char * path ) {
	m_path = std::string( path );
	Volume4dReader();
}

/*
 * destructor
 */
template<class TOutputImage, const unsigned int Dimension>
Volume4dReader<TOutputImage,Dimension>::~Volume4dReader( ) {
}

/*
 * Update
 */
template<class TOutputImage, const unsigned int Dimension>
void
Volume4dReader<TOutputImage,Dimension>::Update( ) {
	typename ArrayReaderType::Pointer reader = ArrayReaderType::New();
	reader->SetFileName( m_path );
	
	typedef itk::CastImageFilter<InputImageType, OutputImageType> 
		CastFilterType;
	typename CastFilterType::Pointer caster = CastFilterType::New();
	caster->SetInput( reader->GetOutput() );
	caster->Update();
	m_output = caster->GetOutput();

}

/*
 * GetOutput
 */
template<class TOutputImage, const unsigned int Dimension>
typename Volume4dReader<TOutputImage,Dimension>::OutputPointerType
Volume4dReader<TOutputImage,Dimension>::GetOutput( ) {
	return m_output;
}

