/**
 * @file BatchMandelCalculator.cc
 * @author David Sladký (xsladk07@stud.fit.vutbr.cz)
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 2022-11-03
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>
#include <cstring>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data = (int *)_mm_malloc(height * width * sizeof(int), ALIGMENT);
	real = (float *)_mm_malloc(BATCH_SIZE * sizeof(float), ALIGMENT);
	img = (float *)_mm_malloc(BATCH_SIZE * sizeof(float), ALIGMENT);
	x = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);
	y = (float *)_mm_malloc(height/2 * sizeof(float), ALIGMENT);
	oneBatchData = (float *)_mm_malloc(BATCH_SIZE * sizeof(float), ALIGMENT);

	initData();
}

BatchMandelCalculator::~BatchMandelCalculator() {
	_mm_free(data);
	_mm_free(real);
	_mm_free(img);
	_mm_free(oneBatchData);
	_mm_free(x);
	_mm_free(y);
	data = NULL;
	real = NULL;
	img = NULL;
	oneBatchData = NULL;
	x = NULL;
	y = NULL;
}


int * BatchMandelCalculator::calculateMandelbrot () {
	
	for(int i = 0; i < height / 2; i++)
	{
		for(int lineBatch = 0; lineBatch < width / BATCH_SIZE; lineBatch++)
		{
			memset(real, 0, BATCH_SIZE * sizeof(float));
			memset(img, 0, BATCH_SIZE * sizeof(float));
			memset(oneBatchData, 0, BATCH_SIZE * sizeof(float));
			float doneCounter = 0;

			for(int l = 0; l < limit && doneCounter < BATCH_SIZE; l++)
			{
				doneCounter = 0;
				#pragma omp simd reduction(+:doneCounter)
				for(int j = 0; j < BATCH_SIZE; j++)
				{
					float r2 = real[j] * real[j];
					float i2 = img[j] * img[j];
					float outOfBound = r2+i2 > 4.0f;

					float cImg = 2.0f * real[j] * img[j] + y[i];
					float cReal = r2 - i2 + x[lineBatch * BATCH_SIZE + j];

					doneCounter += outOfBound;

					oneBatchData[j] = l * ((outOfBound-1.0f)*(-1.0f)) + oneBatchData[j] * outOfBound;

					img[j] = cImg * (outOfBound-1.0f) + img[j] * outOfBound;
					real[j] = cReal * (outOfBound-1.0f) + real[j] * outOfBound;
				}
			}

			int index = i*width + lineBatch * BATCH_SIZE;
			for(int c = 0; c < BATCH_SIZE; c++)
			{
				data[index + c] = static_cast<int>(oneBatchData[c]);
			}
		}

	}

	for(int i = 0; i < height / 2 ; i++)
	{
		memcpy(data+(width)*(height-1)-i*(width), data+i*(width), width * sizeof(int));
	}

	return data;
}

void BatchMandelCalculator::initData()
{
	memset(data, 0, width * height * sizeof(int));
	memset(real, 0, BATCH_SIZE * sizeof(float));
	memset(img, 0, BATCH_SIZE * sizeof(float));
	memset(oneBatchData, 0, BATCH_SIZE * sizeof(float));
	for(int i = 0; i < width; i++)
	{
		x[i] = x_start + i * dx;
	}
	for(int i = 0; i < height/2; i++)
	{
		y[i] = y_start + i * dy;
	}
}