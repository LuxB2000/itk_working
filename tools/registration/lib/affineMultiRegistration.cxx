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
  typename TransformInitializerType::Pointer initializer		= 
                                          TransformInitializerType::New();
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
  metric->SetNumberOfHistogramBins( 50 );
  metric->SetNumberOfSpatialSamples( nbrOfSamples );//5000000 /10 );
  //metric->ReinitializeSeed( 76926294 );
  metric->ReinitializeSeed(  );
	
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

	// configure the initializer
  initializer->SetTransform(   transform );
  initializer->SetFixedImage( m_fixedImage );
  initializer->SetMovingImage( m_movingImage );


	// set the moment instead of the geometry. The gray level moments are
	// better for registration where the anatomical structures aren't centered
	// in the image. The cochlea are not centered in the image.
	//initializer->GeometryOn();
  initializer->MomentsOn();
  initializer->InitializeTransform();

	// init transform with rotation 0
  //transform->SetRotation( rotation );
	itk::Vector<double, 3> initTranslation;
	initTranslation[0] = m_initTranslationX;
	initTranslation[1] = m_initTranslationY;
	initTranslation[2] = m_initTranslationZ;

	transform->Translate( initTranslation );
  registration->SetInitialTransformParameters( transform->GetParameters() );

	// define the stop criteria for the optimizer
	typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
	//number of dimensions: N*(N+1)
	for(unsigned int i=0; i<(Dimension*Dimension); i++ )
	{
		optimizerScales[i] = m_momentScale;
	}
	
	for( unsigned int i=(Dimension*Dimension); 
				i<((Dimension+1)*Dimension); i++ )
	{
		optimizerScales[i] = m_translationScale;
	}

  optimizer->SetScales( optimizerScales );

	// iteration update
  CommandIterationUpdate::Pointer itUpdate = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), itUpdate );
  optimizer->SetNumberOfIterations( m_maxIter );
  optimizer->SetRelaxationFactor( 0.7 );
	optimizer->MinimizeOn();

	// observer
	if( !m_visualize ){
		typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
		typename CommandType::Pointer command = CommandType::New();
		registration->AddObserver( itk::IterationEvent(), command );
	}else{
		/* NEED ITKVTKGLUI
		 * typedef RegistrationVisualizeInterfaceCommand<FixedImageType,
																									MovingImageType,TransformType,
																									RegistrationType,OptimizerType> VisualizeObserverType;
		typename VisualizeObserverType::Pointer command = VisualizeObserverType::New();
		command->SetFixedImage( m_fixedImage );
		command->SetMovingImage( m_movingImage );
		registration->AddObserver( itk::IterationEvent(), command);
		*/
	}

	// set number of level to the pyramid
  registration->SetNumberOfLevels( m_nbrOfLevels );


	//
	// Launch the registration
	//
  try 
    { 
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    } 
  catch( itk::ExceptionObject & err ) 
  {
		throw( err );
  }
  
  m_finalParameters = registration->GetLastTransformParameters();
	m_fixedParameters = transform->GetFixedParameters();

	std::cout << std::endl;
	
}
