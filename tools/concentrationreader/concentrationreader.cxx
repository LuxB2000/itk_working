/*
 * constructors
 */
template<class TOutputImage, const unsigned int Dimension>
ConcentrationReader<TOutputImage,Dimension>::ConcentrationReader( ) {
	
}

template<class TOutputImage, const unsigned int Dimension>
ConcentrationReader<TOutputImage,Dimension>::ConcentrationReader(
																								const char * path ) {
	m_path = std::string( path );
	ConcentrationReader();
}

/*
 * destructor
 */
template<class TOutputImage, const unsigned int Dimension>
ConcentrationReader<TOutputImage,Dimension>::~ConcentrationReader( ) {
}

/*
 * Update
 */
template<class TOutputImage, const unsigned int Dimension>
void
ConcentrationReader<TOutputImage,Dimension>::Update( ) {
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
typename ConcentrationReader<TOutputImage,Dimension>::OutputPointerType
ConcentrationReader<TOutputImage,Dimension>::GetOutput( ) {
	return m_output;
}

