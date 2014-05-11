
#pragma once

#include <stdlib.h>

size_t write_file(char const *filename, unsigned char *file_data, size_t len);
size_t read_file(unsigned char **file_data_ptr, char const *filename);
