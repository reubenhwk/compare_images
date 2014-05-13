
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "compute_asift_keypoints.h"

using namespace std;

void stream_keys_in(ifstream &file, KeyPoints &keys1)
{
	string not_used;
	
	int version;
	file >> not_used >> version;
	// Follow the same convention of David Lowe:
	// the first line contains the number of keypoints and the length of the desciptors (128)
	int len;
	file >> not_used >> len;

	int tilts = 0;

	file >> not_used >> keys1.tilts;
	keys1.keys.resize(keys1.tilts);
	for (int tt = 0; tt < keys1.tilts; tt++) {
		int rotations = 0;
		file >> not_used >> rotations;
		keys1.keys[tt].resize(rotations);
		for (int rr = 0; rr < rotations; rr++) {
			int keys = 0;
			file >> not_used >> keys;
			for (int i = 0; i < keys; i++) {
				keypoint kp;
				file >> kp.x;
				file >> kp.y;
				file >> kp.scale;
				file >> kp.angle;

				for (int ii = 0; ii < len; ii++) {
					file >> kp.vec[ii];
				}
				keys1.keys[tt][rr].push_back(kp);
				++keys1.count;
			}
		}
	}
}

void stream_keys_out(ofstream &file, KeyPoints &keys1)
{
	/* set width to 3 to reserve the first three bytes
	 * for the version number.  The fourth byte is always
	 * for the newline. */
	file << std::setw(10) << "version:" << ' ' << std::setw(10) << 1 << std::endl;

	file << std::setw(10) << "length:" << ' ' << std::setw(10) << VecLength << std::endl;

	file << std::setw(10) << "tilts:" << ' ' << std::setw(10) << (int)keys1.keys.size() << std::endl;
	for (int tt = 0; tt < (int)keys1.keys.size(); tt++) {
		file << std::setw(10) << "rotations:" << ' ' << std::setw(10) << (int)keys1.keys[tt].size() << std::endl;
		for (int rr = 0; rr < (int)keys1.keys[tt].size(); rr++) {
			keypointslist::iterator ptr = keys1.keys[tt][rr].begin();
			file << std::setw(10) << "keys:" << ' ' << std::setw(10) << (int)keys1.keys[tt][rr].size() << std::endl;
			for (int i = 0; i < (int)keys1.keys[tt][rr].size(); i++, ptr++) {
				file << std::setw(10) << ptr->x << ' ';
				file << std::setw(10) << ptr->y << ' ';
				file << std::setw(10) << ptr->scale << ' ';
				file << std::setw(10) << ptr->angle;

				for (int ii = 0; ii < (int)VecLength; ii++) {
					file << ' ' << std::setw(10) << ptr->vec[ii];
				}

				file << std::endl;
			}
		}
	}
}

