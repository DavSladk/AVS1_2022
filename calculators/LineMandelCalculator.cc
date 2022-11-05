/**
 * @file LineMandelCalculator.cc
 * @author David Sladký (xsladk07@stud.fit.vutbr.cz)
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 2022-11-03
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
	real = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);
	img = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);
	x = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);
	y = (float *)_mm_malloc(height/2 * sizeof(float), ALIGMENT);
	oneLineData = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);

	initData();
}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	_mm_free(real);
	_mm_free(img);
	_mm_free(oneLineData);
	_mm_free(x);
	_mm_free(y);
	data = NULL;
	real = NULL;
	img = NULL;
	oneLineData = NULL;
	x = NULL;
	y = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	for(int i = 0; i < height / 2 ; i++)
	{
		memset(real, 0, width * sizeof(float));
		memset(img, 0, width * sizeof(float));
		memset(oneLineData, 0, width * sizeof(float));
		float doneCounter = 0;

		for(int l = 0; l < limit && doneCounter < width; l++)
		{
			doneCounter = 0.0f;
			#pragma omp simd reduction(+:doneCounter)
			for(int j = 0; j < width; j++)
			{
				float r2 = real[j] * real[j];
				float i2 = img[j] * img[j];
				float outOfBound = r2+i2 > 4.0f;

				float cImg = 2.0f * real[j] * img[j] + y[i];
				float cReal = r2 - i2 + x[j];

				doneCounter += outOfBound;

				oneLineData[j] = l * ((outOfBound-1.0f)*(-1.0f)) + oneLineData[j] * outOfBound;

				img[j] = cImg * (outOfBound-1.0f) + img[j] * outOfBound;
				real[j] = cReal * (outOfBound-1.0f) + real[j] * outOfBound;
			}
		}
		
		int index = i*width;
		for(int c = 0; c < width; c++)
		{
			data[index + c] = static_cast<int>(oneLineData[c]);
		}
		memcpy(data+(width)*(height-1)-i*(width), data+i*(width), width * sizeof(int));
	}

	return data;
}

void LineMandelCalculator::initData()
{
	memset(data, 0, width * height * sizeof(int));
	memset(real, 0, width * sizeof(float));
	memset(img, 0, width * sizeof(float));
	memset(oneLineData, 0, width * sizeof(float));
	for(int i = 0; i < width; i++)
	{
		x[i] = x_start + i * dx;
	}
	for(int i = 0; i < height/2; i++)
	{
		y[i] = y_start + i * dy;
	}
}