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
	real = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);
	img = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);
	x = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);
	y = (float *)_mm_malloc(height * sizeof(float), ALIGMENT);
	isDone = (bool *)_mm_malloc(width * sizeof(bool), ALIGMENT);
	
	initData();
}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	_mm_free(real);
	_mm_free(img);
	_mm_free(isDone);
	_mm_free(x);
	_mm_free(y);
	data = NULL;
	real = NULL;
	img = NULL;
	isDone = NULL;
	x = NULL;
	y = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	for(int i = 0; i < height / 2 ; i++)
	{
		memset(real, 0, width * sizeof(float));
		memset(img, 0, width * sizeof(float));
		memset(isDone, 0, width * sizeof(bool));
		int doneCounter = 0;
		
		for(int l = 0; l < limit && doneCounter < width; l++)
		{
			doneCounter = 0;
			#pragma omp simd reduction(+:doneCounter)
			for(int j = 0; j < width; j++)
			{
				int index = i * width + j;

				float r2 = real[j] * real[j];
				float i2 = img[j] * img[j];
				bool outOfBound = r2+i2 > 4;

				float cImg = 2.0f * real[j] * img[j] + y[i];
				float cReal = r2 - i2 + x[j];

				isDone[j] = isDone[j] || outOfBound;
				doneCounter += isDone[j];

				data[index] = l * (!isDone[j]) + data[index] * isDone[j];

				img[j] = cImg * (!isDone[j]) + img[j] * isDone[j];
				real[j] = cReal * (!isDone[j]) + real[j] * isDone[j];
			}
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
	memset(isDone, 0, width * sizeof(bool));
	for(int i = 0; i < width; i++)
	{
		x[i] = x_start + i * dx;
	}
	for(int i = 0; i < height; i++)
	{
		y[i] = y_start + i * dy;
	}
}