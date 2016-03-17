/*
 * =====================================================================================
 *
 *       Filename:  writeIndexFile.h
 *
 *    Description:  Write index into a specified file, Note: it writes at the end of the 
 *    file if it is already exist.
 *
 *        Version:  1.0
 *        Created:  30/04/15 11:55:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

#ifndef __WRITEINDEXFILE__H
#define __WRITEINDEXFILE__H

#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "itkMacro.h"

template <class TIndexType, const unsigned int Dimension=3>
class WriteIndexFile{
public:
	// constructor and destructor
	WriteIndexFile(){
		m_separator = char(',');
	}
	~WriteIndexFile(){
		//if( m_indexList ) delete m_indexList;
	}
	

	typedef TIndexType IndexType;
	typedef std::vector< IndexType > VectorIndexType;

	// get the index list 
	void SetIndexList(VectorIndexType v){
		m_indexList = VectorIndexType(v);
	}

	// set the separator in the index file, default = ','
	void SetSeparator( char s ){
		m_separator = s;
	}

	// set the file to read
	void SetFileName(std::string fName ){
		m_filePath = fName;
	}

	// use this function to write in the file
	void Update(){
		if( m_filePath.length() == 0 ){
			itk::ExceptionObject err("Set the path to reach the file.");
			throw( err );
		}

		std::ofstream myfile ( m_filePath.c_str(), std::ofstream::app );

		if( !myfile.is_open() ){
			myfile.close();
			itk::ExceptionObject err("Error while opening the file: " + m_filePath  );
			throw( err );
		}

		std::string num;
		unsigned int i=0, d=0, nbrOfIndex = m_indexList.size();
		IndexType ind;

		while( i<nbrOfIndex ){
			ind = m_indexList.at(i);
			d = 0;
			while( d<Dimension ){
				myfile << ind[d] ;
				//if( !( i==(nbrOfIndex-1) & (d==Dimension-1) ) ){
					myfile << m_separator;
				//}
				d++;
			}
			i++;
		}

		//myfile << "\n";
		myfile.close();
	} // end update

private:
	VectorIndexType m_indexList;
	std::string m_filePath;
	char m_separator;
};


#endif

