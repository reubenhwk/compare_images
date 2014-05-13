// Copyright (c) 2008-2011, Guoshen Yu <yu@cmap.polytechnique.fr>
// Copyright (c) 2008-2011, Jean-Michel Morel <morel@cmla.ens-cachan.fr>
//
// WARNING:
// This file implements an algorithm possibly linked to the patent
//
// Jean-Michel Morel and Guoshen Yu, Method and device for the invariant
// affine recognition recognition of shapes (WO/2009/150361), patent pending.
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.
//
//
//*----------------------------- demo_ASIFT  --------------------------------*/
// Detect corresponding points in two images with the ASIFT method.

// Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
//
// Reference: J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image
//            Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009.
// Reference: ASIFT online demo (You can try ASIFT with your own images online.)
//                        http://www.ipol.im/pub/algo/my_affine_sift/
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "demo_lib_sift.h"
#include "lib/io_png/io_png.h"
#include "lib/magickio/magickio.h"

#include "library.h"
#include "frot.h"
#include "fproj.h"
#include "compute_asift_keypoints.h"
#include "compute_asift_matches.h"

#include "image_size.h"
#include "keyio.h"

int main(int argc, char **argv)
{

	if (argc != 3) {
		std::cerr << " ******************************************************************************* " << std::endl
		    << " ***************************  ASIFT image matching  **************************** " << std::endl
		    << " ******************************************************************************* " << std::endl
		    << "Usage: " << argv[0] << " key1.txt key2.txt" << std::endl
		    << " ******************************************************************************* " << std::endl
		    << " *********************  Jean-Michel Morel, Guoshen Yu, 2010 ******************** " << std::endl
		    << " ******************************************************************************* " << std::endl;
		return 1;
	}

	KeyPoints keys[2];
	for (int i = 0; i < 2; ++i) {
		std::ifstream file_key(argv[i + 1]);
		if (file_key.is_open()) {
			stream_keys_in(file_key, keys[i]);
		} else {
			std::cerr << "Unable to open the file keys2 for reading." << std::endl;
		}
		file_key.close();
	}

	//// Match ASIFT keypoints
	matchingslist matchings;
	cout << "Matching the keypoints..." << endl;
	time_t const tstart = time(0);

	// Define the SIFT parameters
	siftPar siftparameters;
	default_sift_parameters(siftparameters);

	float wS1, hS1, wS2, hS2;
	wS1 = hS1 = wS2 = hS2 = 1;
	
	int rc = compute_asift_matches(keys[0].tilts, keys[1].tilts, 
						keys[0].width, keys[0].height,
						keys[1].width, keys[1].height,
						0,
						keys[0].keys, keys[1].keys,
						matchings, siftparameters);
	time_t const tend = time(0);
	cout << "Keypoints matching accomplished in " << difftime(tend, tstart) << " seconds." << endl;

	return rc;
}
