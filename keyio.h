
#pragma once

#include <iostream>
#include <vector>

#include "compute_asift_keypoints.h"

void stream_keys_in(std::ifstream &file1, int &num_keys1, double &zoom1, vector < vector < keypointslist > > &keys1);
void stream_keys_out(std::ofstream &file1, int num_keys1, double zoom1, vector < vector < keypointslist > > &keys1);

