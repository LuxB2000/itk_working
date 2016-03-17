/*
 * =====================================================================================
 *
 *       Filename:  debug.h
 *
 *    Description:  debug function
 *
 *        Version:  1.0
 *        Created:  23/04/15 13:37:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

#ifndef __DEBUG__H
#define __DEBUG__H

#include "color.h"

#define DEBUG 1

void debug( std::string msg ){
	if( DEBUG ){
		blueMessage( msg );
	}
}

#endif
