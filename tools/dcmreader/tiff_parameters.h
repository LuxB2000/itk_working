#ifndef __TIFFPARAMETERS__
#define __TIFFPARAMETERS__

/*
 * Contains the tiff parameters.
 *
 */

#include <vector>
#include <string>
#include <time.h>

class TiffParameters
{
	public:
		typedef std::vector<double> VectorType;
		TiffParameters()
		{
			spacing = VectorType(3,1); // size 3 with default values 1
			origin = VectorType(3,0);
			orientation = VectorType(3,0);
		}
		~TiffParameters()
		{
		}

		VectorType spacing;
		VectorType origin;
		VectorType orientation;
		double tr, te, ti, fliplist;
		struct tm time_run, time_complete;

	private:
};


#endif
