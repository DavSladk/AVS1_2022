/**
 * @file LineMandelCalculator.cc
 * @author David Sladký (xsladk07@stud.fit.vutbr.cz)
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 2022-10-30
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <cstring>


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	data = (int *)_mm_malloc(height * width * sizeof(int), ALIGMENT);
	real = (float *)_mm_malloc(height * width * sizeof(int), ALIGMENT);
	img = (float *)_mm_malloc(height * width * sizeof(int), ALIGMENT);
	initData();
}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	_mm_free(real);
	_mm_free(img);
	data = NULL;
	real = NULL;
	img = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	for(int i = 0; i < height / 2 ; i++)
	{
		float y = y_start + i * dy; // current imaginary value
		int doneCount = 0;
		for(int l = 0; l < limit; l++)
		{
			#pragma omp simd reduction(+:doneCount)//aligned(data, real, img: 64)
			for(int j = 0; j < width; j++)
			{
				float x = x_start + j * dx; // current real value

				int index = i * width + j;

				float r2 = real[index] * real[index];
				float i2 = img[index] * img[index];

				if(r2+i2 > 4)
				{
					data[index] = std::min(data[index], l);
				}
				else
				{
					img[index] = 2.0f * real[index] * img[index] + y;
					real[index] = r2 - i2 + x;
				}
			}
		}
		//memcpy(data+(width)*(height-1)-i*(width), data+i*(width), width * sizeof(int));
	}
	// memcpy is outside of main loop, because compiler doesn't like it inside
	for(int i = 0; i < height / 2 ; i++)
	{
		memcpy(data+(width)*(height-1)-i*(width), data+i*(width), width * sizeof(int));
	}

	return data;
}

void LineMandelCalculator::initData()
{
	for(int i = 0; i < height * width; i++)
	{
		data[i] = limit;
	}
	memset(real, 0, height * width * sizeof(float));
	memset(img, 0, height * width * sizeof(float));

}