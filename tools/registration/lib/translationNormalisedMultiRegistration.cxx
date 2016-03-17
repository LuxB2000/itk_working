/*
 * Constructor
 */
template<class TMovingImage, class TFixedImage, const unsigned int Dimension>
AffineMultiRegistration<TMovingImage,TFixedImage,Dimension>::AffineMultiRegistration(){

	// default values
	m_pixelRatio = 40;
	m_nbrOfLevels = 3;
	m_momentScale = 1.0;
	m_translationScale = (double) 1/1E5;
	m_initTranslationX = 0;
	m_initTranslationY = 0;
	m_initTranslationZ = 0;
	m_maxIter = 200;
	m_visualize = false; //TODO Must be false

}
/*
 * Constructor
 */
template<class TMovingImage, class TFixedImage, const unsigned int Dimension>
AffineMultiRegistration<TMovingImage,TFixedImage,Dimension>::~AffineMultiRegistration(){
}

/*
 * Update
 */
template<class TMovingImage, class TFixedImage, const unsigned int Dimension>
void
AffineMultiRegistration<TMovingImage,TFixedImage,Dimension>::Update(){

	//
	// instantiation
	//
  typename MetricType::Pointer         metric        = MetricType::New();
  typename OptimizerType::Pointer      optimizer     = OptimizerType::New();
  typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  typename RegistrationType::Pointer   registration  = RegistrationType::New();
  typename TransformType::Pointer			transform			= TransformType::New();
  typename FixedImagePyramidType::Pointer fixedImagePyramid = 
																					FixedImagePyramidType::New();
  typename MovingImagePyramidType::Pointer movingImagePyramid =
																					MovingImagePyramidType::New();
	
	//
	// metric
	//
	const long nbrOfSamples = (long) ( 
		m_fixedImage->GetLargestPossibleRegion().GetSize()[0] *
		m_fixedImage->GetLargestPossibleRegion().GetSize()[1] *
		m_fixedImage->GetLargestPossibleRegion().GetSize()[2] *
		m_pixelRatio / 100 );
	std::cout<<"Registration is done with "<< nbrOfSamples 
		<<" ("<<m_pixelRatio<< " percent of pixels) and "
		<<m_nbrOfLevels<<" levels."<<std::endl;
  //metric->SetNumberOfHistogramBins( 248 );
  metric->SetNumberOfSpatialSamples( nbrOfSamples );//5000000 /10 );
  //metric->ReinitializeSeed( 76926294 );
  //metric->ReinitializeSeed(  );
	metric->SetFixedImageStandardDeviation( 0.4 );
	metric->SetMovingImageStandardDeviation( 0.4 );
	
	//
	// pipeline connections
	//
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  registration->SetTransform(     transform     );
  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );


	// connect the inputs
  registration->SetFixedImage( m_fixedImage   );
  registration->SetMovingImage( m_movingImage );

	// set the regions
  registration->SetFixedImageRegion( 
     m_fixedImage->GetLargestPossibleRegion() );


	// iteration update
  CommandIterationUpdate::Pointer itUpdate = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), itUpdate );
  optimizer->SetNumberOfIterations( m_maxIter );
  //optimizer->SetRelaxationFactor( 0.7 );
	optimizer->SetLearningRate( 1 );
	optimizer->MaximizeOn();

		typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
		typename CommandType::Pointer command = CommandType::New();
		registration->AddObserver( itk::IterationEvent(), command );

	// set number of level to the pyramid
  registration->SetNumberOfLevels( m_nbrOfLevels );

	typedef typename RegistrationType::ParametersType ParametersType;
	ParametersType initParam( transform->GetNumberOfParameters() );
	// rotation matrix
	initParam.Fill( 0.0 ); // init translation and rotation matrix

	// init translation
	initParam[0] = m_initTranslationX;
	initParam[1] = m_initTranslationY;
	initParam[2] = m_initTranslationZ;
	
	registration->SetInitialTransformParameters( initParam );
	
	//
	// Launch the registration
	//
  try 
    { 
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
		m_finalParameters = registration->GetLastTransformParameters();
		m_fixedParameters = transform->GetFixedParameters();
    } 
  catch( itk::ExceptionObject & err ) 
  {
		throw( err );
  }
  

	std::cout << std::endl;
	
}
