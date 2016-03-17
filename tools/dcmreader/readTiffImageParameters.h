#ifndef __READTIFFIMAGEPARAMETERS__
#define __READTIFFIMAGEPARAMETERS__

/*
 *
 *
 * The file is assume to be formated as followed:
 * -- This is a comment
 * <NAME> = <VALUE>
 */
#include "tiff_parameters.h"
#include <fstream>
#include <string>

class ReadTiffImageParameters
{
	public:
		ReadTiffImageParameters(std::string path);
		~ReadTiffImageParameters();
		void Read();
		TiffParameters* GetParamters(){ return m_imageParameters; }
		typedef std::vector<int> DateVectorType;

	private:
		std::string m_path;
		TiffParameters* m_imageParameters;
		void fromDateStrToDateStruct(std::string, std::vector<int>*);
};


#ifndef ITK_MANUAL_INSTANTIATION
#include "readTiffImageParameters.cxx"
#endif

#endif
