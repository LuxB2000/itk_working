#include "readTiffImageParameters.h"

ReadTiffImageParameters::ReadTiffImageParameters( std::string path )
{
	m_path = path;
	m_imageParameters = new TiffParameters();
}

ReadTiffImageParameters::~ReadTiffImageParameters()
{
	if(m_imageParameters)
		delete m_imageParameters;
}

void ReadTiffImageParameters::fromDateStrToDateStruct( std::string dateStr, 
		ReadTiffImageParameters::DateVectorType* vect ){
	//std::cout << "TEST : " << dateStr << std::endl;
	// Ex: "20100628T150721"
	vect->at(0) = atoi(dateStr.substr(1,4).c_str()); // year
	vect->at(2) = atoi(dateStr.substr(5,2).c_str()); // day
	vect->at(1) = atoi(dateStr.substr(7,2).c_str()); // month
	vect->at(3) = atoi(dateStr.substr(10,2).c_str()); // hour
	vect->at(4) = atoi(dateStr.substr(12,2).c_str()); // min
	vect->at(5) = atoi(dateStr.substr(14,2).c_str()); // sec

	//std::cout << "TEST2: " <<
	//	"y= " << vect->at(0) << " " << dateStr.substr(1,4) << "\n" <<
	//	"m= " << vect->at(2) << " " << dateStr.substr(5,2) << "\n" <<
	//	"d= " << vect->at(1) << " "<< dateStr.substr(7,2) << "\n" <<
	//	"h= " << vect->at(3) << " "<< dateStr.substr(10,2) << "\n" <<
	//	"m= " << vect->at(4) << " "<< dateStr.substr(12,2) << "\n" <<
	//	"s= " << vect->at(5) << " " << dateStr.substr(14,2) << "\n" <<
	//	std::endl;


}

void
ReadTiffImageParameters::Read()
{
	double lro=1, np=2, lpe=1, nv=1, thk=1, lpe2=0,
				 nv2=0, pro=0, ppe=0, pss0=0, tr_=0.0,
				 te=0.0, ti=0.0, fliplist=0;
	std::ifstream file(m_path.c_str());
	DateVectorType tr, tc;

	// parse the file
	for( std::string line; getline(file, line); )
	{
		// check if it is a comment
		if( line.size() > 1 ) // empty line
		if( line.substr(0,2).compare(std::string("--")) != 0 )
		{
			// it's not a comment
			// get the name of the parameters <NAME> = <VALUE>
			std::string name = line.substr(0,line.find_first_of("=")-1);
			//std::cout<<"\t<NAME>="<<name;

			// get the value
			std::string value = line.substr(line.find_first_of("=")+2);
			//remove potential space in value
			value.erase(std::remove(value.begin(),value.end(),' '),value.end());
				
			//std::cout<<" <VALUE>="<<value<<std::endl;

			if( name.compare("date") == 0 )
			{
			}else if( name.compare("time_run") == 0 ){
				tr = DateVectorType(6,0);
				fromDateStrToDateStruct( value, &tr );
			}else if( name.compare("time_complete") == 0 ){
				tc = DateVectorType(6,0);
				fromDateStrToDateStruct( value, &tc );
			}else if( name.compare("lro") == 0 ){
				lro=atof(value.c_str());
			}else if( (name.compare("lpe") == 0) || (name.compare("lpe0") == 0 )){
				lpe=atof(value.c_str());
			}else if( name.compare("lpe2") == 0 ){
				lpe2=atof(value.c_str());
			}else if( name.compare("np") == 0 ){
				np=atof(value.c_str());
			}else if( (name.compare("nv") == 0) || (name.compare("nv0") == 0) ){
				nv=atof(value.c_str());
			}else if( name.compare("nv2") == 0 ){
				nv2=atof(value.c_str());
			}else if( name.compare("thk") == 0 ){
				thk = atof( value.c_str() );
			}else if( name.compare("pro") == 0 ){
				pro=atof(value.c_str());
			}else if( name.compare("ppe") == 0 ){
				ppe=atof(value.c_str());
			}else if( name.compare("pss0") == 0 ){
				pss0=atof(value.c_str());
			}else if( name.compare("tr") == 0 ){
				tr_=atof(value.c_str());
			}else if( name.compare("ti") == 0 ){
				ti=atof(value.c_str());
			}else if( name.compare("te") == 0 ){
				te=atof(value.c_str());
			}else if( name.compare("fliplist") == 0 ){
				fliplist = atof(value.c_str());
			}
			else //default:
			{
				// REMOVE IT FOR DEBUGING
				if(0){
				std::cerr<<"Warning parameters with Name <"
					<<name<<"> is not recognised.";
				 std::cerr<<" Value is <"<<value<<">.\n";
				std::cerr<<"We ignore it.\n"<<std::endl;
				}
			}

		}
	}

	// set parameters
	m_imageParameters->spacing.at(0) = lpe/nv * 10; // factor 10 to reach mm
	m_imageParameters->spacing.at(1) = lro/(np*0.5) * 10;
	if( (lpe2 != 0) && (nv2 != 0) )
	{
		m_imageParameters->spacing.at(2) = lpe2/nv2;
	}else{
		m_imageParameters->spacing.at(2) = thk;
	}

	// std::cout<<"["<<pro<<","<<ppe<<","<<pss0<<"]"<<std::endl;
	m_imageParameters->origin.at(0) = pro;
	m_imageParameters->origin.at(1) = ppe;
	m_imageParameters->origin.at(2) = pss0;
	
	// set the time parameters
	m_imageParameters->time_run.tm_year = tr.at(0);
	m_imageParameters->time_run.tm_mon = tr.at(1);
	m_imageParameters->time_run.tm_mday = tr.at(2);
	m_imageParameters->time_run.tm_hour = tr.at(3);
	m_imageParameters->time_run.tm_min = tr.at(4);
	m_imageParameters->time_run.tm_sec = tr.at(5);

	m_imageParameters->time_complete.tm_year = tc.at(0);
	m_imageParameters->time_complete.tm_mon = tc.at(1);
	m_imageParameters->time_complete.tm_mday = tc.at(2);
	m_imageParameters->time_complete.tm_hour = tc.at(3);
	m_imageParameters->time_complete.tm_min = tc.at(4);
	m_imageParameters->time_complete.tm_sec = tc.at(5);

	m_imageParameters->tr = tr_;
	m_imageParameters->te = te;
	m_imageParameters->ti = ti;
	m_imageParameters->fliplist = fliplist;

			//std::cout<<"Test: "<<m_imageParameters->spacing.at(0)<<"  "<<m_imageParameters->spacing.at(1)<<"  "<<m_imageParameters->spacing.at(2)<<std::endl;
}
