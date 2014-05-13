
#include <stdio.h>
#include <stdlib.h>

size_t read_file(unsigned char **file_data_ptr, char const *filename)
{
	FILE *in = fopen(filename, "rb");

	if (!in) {
		return 0;
	}

	fseek(in, 0, SEEK_END);
	size_t allocated = ftell(in);
	unsigned char *file_data = malloc(allocated);

	rewind(in);
	size_t used = fread(file_data, 1, allocated, in);
	while (used > 0 && !feof(in)) {
		if (used >= allocated) {
			allocated += 4096;
			file_data = realloc(file_data, allocated);
		}
		used += fread(file_data + used, 1, allocated - used, in);
	}

	fclose(in);
	*file_data_ptr = file_data;

	return used;
}

size_t write_file(char const *filename, unsigned char *file_data, size_t len)
{
	FILE *out = fopen(filename, "wb");

	if (!out) {
		return (size_t) - 1;
	}

	size_t count = fwrite(file_data, 1, len, out);

	fclose(out);

	return count;
}
