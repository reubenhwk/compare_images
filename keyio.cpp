
#include <fstream>
#include <iomanip>
#include <vector>

#include "compute_asift_keypoints.h"

using namespace std;

void stream_keys_in(ifstream &file1, int &num_keys1, double &zoom1, vector < vector < keypointslist > > &keys1)
{
	int version;
	file1 >> version;
	// Follow the same convention of David Lowe:
	// the first line contains the number of keypoints and the length of the desciptors (128)
	file1 >> num_keys1;
	int len;
	file1 >> len;

	for (int tt = 0; tt < (int)keys1.size(); tt++) {
		for (int rr = 0; rr < (int)keys1[tt].size(); rr++) {
			keypointslist::iterator ptr = keys1[tt][rr].begin();
			for (int i = 0; i < (int)keys1[tt][rr].size(); i++, ptr++) {
				file1 >> ptr->x;
				file1 >> ptr->y;
				file1 >> ptr->scale;
				file1 >> ptr->angle;

				for (int ii = 0; ii < len; ii++) {
					file1 >> ptr->vec[ii];
				}
			}
		}
	}
}

void stream_keys_out(ofstream &file1, int num_keys1, double zoom1, vector < vector < keypointslist > > &keys1)
{
	/* set width to 3 to reserve the first three bytes
	 * for the version number.  The fourth byte is always
	 * for the newline. */
	file1 << std::setw(10) << 1 << std::endl;
	// Follow the same convention of David Lowe:
	// the first line contains the number of keypoints
	file1 << std::setw(10) << num_keys1 << std::endl;

	file1 << std::setw(10) << "keys1.size() " << (int)keys1.size() << std::endl;
	for (int tt = 0; tt < (int)keys1.size(); tt++) {
		file1 << "keys1[tt].size() " << std::setw(10) << (int)keys1[tt].size() << std::endl;
		for (int rr = 0; rr < (int)keys1[tt].size(); rr++) {
			keypointslist::iterator ptr = keys1[tt][rr].begin();
			file1 << "keys1[tt][rr].size() " << std::setw(10) << (int)keys1[tt][rr].size() << std::endl;
			for (int i = 0; i < (int)keys1[tt][rr].size(); i++, ptr++) {
				file1 << std::setw(10) << zoom1 * ptr->x << ' ';
				file1 << std::setw(10) << zoom1 * ptr->y << ' ';
				file1 << std::setw(10) << zoom1 * ptr->scale << ' ';
				file1 << std::setw(10) << ptr->angle;

				file1 << std::setw(10) << VecLength;
				for (int ii = 0; ii < (int)VecLength; ii++) {
					file1 << ' ' << std::setw(10) << ptr->vec[ii];
				}

				file1 << std::endl;
			}
		}
	}
}


