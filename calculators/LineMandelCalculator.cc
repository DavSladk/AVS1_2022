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

//#define REDUCE

LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	// @TODO allocate & prefill memory
	data = (int *)(aligned_alloc(ALIGMENT, height * width * sizeof(int)));
	real = (float *)(aligned_alloc(ALIGMENT, height * width * sizeof(float)));
	img = (float *)(aligned_alloc(ALIGMENT, height * width * sizeof(float)));
	initData();

}

LineMandelCalculator::~LineMandelCalculator() {
	// @TODO cleanup the memory
	free(data);
	free(real);
	free(img);
	data = NULL;
	real = NULL;
	img = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	// @TODO implement the calculator & return array of integers
	int* pdata = data;
	float* preal = real;
	float* pimg = img;
	
	
	for(int i = 0; i < height / 2 + 1; i++)
	{
		#ifndef REDUCE
		float y = y_start + i * dy; // current imaginary value

		for(int l = 0; l < limit; l++)
		{
			#pragma omp simd aligned(pdata, preal, pimg: 64)
			for(int j = 0; j < width; j++)
			{
				float x = x_start + j * dx; // current real value

				int currentIndex = i * width + j;
				float r2 = preal[currentIndex] * preal[currentIndex];
				float i2 = pimg[currentIndex] * pimg[currentIndex];
				
				if(r2 + i2 > 4.0f && pdata[currentIndex] == limit)
				{
					pdata[currentIndex] = l;
				}
				else if(r2 + i2 < 4.0f)
				{
					pimg[currentIndex] = 2.0f * preal[currentIndex] * pimg[currentIndex] + y;
					preal[currentIndex] = r2 - i2 + x;
				}

			}
		}
		#else
		int doneCount = 0;
		float y = y_start + i * dy; // current imaginary value

		for(int l = 0; (l < limit) && (doneCount < limit); l++)
		{
			#pragma omp simd aligned(pdata, preal, pimg: 64) reduction(+:doneCount)
			for(int j = 0; j < width; j++)
			{
				float x = x_start + j * dx; // current real value

				int currentIndex = i * width + j;
				float r2 = preal[currentIndex] * preal[currentIndex];
				float i2 = pimg[currentIndex] * pimg[currentIndex];
				
				if(r2 + i2 > 4.0f && pdata[currentIndex] == limit)
				{
					pdata[currentIndex] = l;
					doneCount += 1;
				}
				else if(r2 + i2 < 4.0f)
				{
					pimg[currentIndex] = 2.0f * preal[currentIndex] * pimg[currentIndex] + y;
					preal[currentIndex] = r2 - i2 + x;
					doneCount += 0;
				}

			}
		}
		#endif

		memcpy(data+(width)*(height)-i*(width), data+i*(width), width * sizeof(int));
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