// Authors: Unknown. Please, if you are the author of this file, or if you
// know who are the authors of this file, let us know, so we can give the
// adequate credits and/or get the adequate authorizations.

#include "flimage.h"

#include <string.h>

//////////////////////////////////////////////// Class flimage
//// Construction
flimage::flimage():width(0), height(0), p(0)
{
}

flimage::flimage(int w, int h):width(w), height(h), p(new float[w * h])
{
	for (int j = width * height - 1; j >= 0; j--)
		p[j] = 0.0;
}

flimage::flimage(int w, int h, float v):width(w), height(h), p(new float[w * h])
{
	for (int j = width * height - 1; j >= 0; j--)
		p[j] = v;
}

flimage::flimage(int w, int h, float const *v):width(w), height(h), p(new float[w * h])
{
	memcpy(p, v, sizeof(float) * w * h);
}

void flimage::create(int w, int h)
{
	erase();
	width = w;
	height = h;
	p = new float[w * h];
	memset(p, 0, sizeof(float) * w * h);
}

void flimage::create(int w, int h, float const *v)
{
	erase();
	width = w;
	height = h;
	p = new float[w * h];
	memcpy(p, v, sizeof(float) * w * h);
}

flimage::flimage(const flimage & im):width(im.width), height(im.height), p(new float[im.width * im.height])
{
	memcpy(p, im.p, sizeof(float) * width * height);
}

flimage & flimage::operator=(const flimage & im)
{
	if (&im == this) {
		return *this;
	}

	if (width != im.width || height != im.height) {
		erase();
		width = im.width;
		height = im.height;
		p = new float[width * height];
	}

	memcpy(p, im.p, sizeof(float) * width * height);

	return *this;
}

//// Destruction
void flimage::erase()
{
	width = height = 0;
	if (p)
		delete[]p;
	p = 0;
}

flimage::~flimage()
{
	erase();
}
