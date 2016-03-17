/*
 * =====================================================================================
 *
 *       Filename:  readIndexFile.h
 *
 *    Description:  A short class that read list of index and return a std::vector.
 *		Default usage: separator is ','. Use SetSeparator for a change
 *
 *        Version:  1.0
 *        Created:  30/04/15 11:02:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

#ifndef __READINDEXFILE__H
#define __READINDEXFILE__H

#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "itkMacro.h"

template <class TIndexType, const unsigned int Dimension=3>
class ReadIndexFile{
public:
	// constructor and destructor
	ReadIndexFile(){
		m_separator = char(',');
	}
	~ReadIndexFile(){
		//if( m_indexList ) delete m_indexList;
	}
	

	typedef TIndexType IndexType;
	typedef std::vector< IndexType > VectorIndexType;

	// get the index list 
	VectorIndexType GetIndexList(){
		return m_indexList;
	}

	// set the separator in the index file, default = ','
	void SetSeparator( char s ){
		m_separator = s;
	}

	// set the file to read
	void SetFileName(std::string fName ){
		m_filePath = fName;
	}

	// use this function to read the file
	void Update(){
		if( m_filePath.length() == 0 ){
			itk::ExceptionObject err("Set the path to reach the file.");
			throw( err );
		}

		std::ifstream myfile ( m_filePath.c_str() );

		if( !myfile.good() ){
			myfile.close();
			itk::ExceptionObject err("Error while reading the file: " + m_filePath  );
			throw( err );
		}

		std::string num;
		unsigned int i=0;
		IndexType ind;
		m_indexList = VectorIndexType();

		while( std::getline( myfile, num, m_separator) ){
			//std::cout << num << std::endl;
			ind[i] = atof( num.c_str() );
			i++;
			if( i==(Dimension) ){
				i=0;
				m_indexList.push_back( ind );
				//std::cout << ind << std::endl;
			}
		}


		myfile.close();
	} // end update

private:
	VectorIndexType m_indexList;
	std::string m_filePath;
	char m_separator;
};


#endif
