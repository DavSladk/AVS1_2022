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
	real = (float *)_mm_malloc(BLOCK_SIZE * sizeof(float), ALIGMENT);
	img = (float *)_mm_malloc(BLOCK_SIZE * sizeof(float), ALIGMENT);
	x = (float *)_mm_malloc(width * sizeof(float), ALIGMENT);
	y = (float *)_mm_malloc(height * sizeof(float), ALIGMENT);
	isDone = (bool *)_mm_malloc(BLOCK_SIZE * sizeof(bool), ALIGMENT);
	
	initData();
}

BatchMandelCalculator::~BatchMandelCalculator() {
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


int * BatchMandelCalculator::calculateMandelbrot () {
	for(int blockY = 0; blockY < height / (2 * BLOCK_SIZE); blockY++)
	{
		int blockYindex = blockY * BLOCK_SIZE * width;
		for(int blockX = 0; blockX < width / BLOCK_SIZE; blockX++)
		{
			int blockXindex = blockX * BLOCK_SIZE;
			int blockXYindex = blockYindex + blockXindex;
			for(int i = 0; i < BLOCK_SIZE; i++)
			{
				memset(real, 0, BLOCK_SIZE * sizeof(float));
				memset(img, 0, BLOCK_SIZE * sizeof(float));
				memset(isDone, 0, BLOCK_SIZE * sizeof(bool));
				int doneCounter = 0;

				for(int l = 0; l < limit && doneCounter < BLOCK_SIZE; l++)
				{
					doneCounter = 0;
					#pragma omp simd reduction(+:doneCounter)
					for(int j = 0; j < BLOCK_SIZE; j++)
					{
						int index = blockXYindex + j + i*width;

						float r2 = real[j] * real[j];
						float i2 = img[j] * img[j];
						bool outOfBound = r2+i2 > 4;

						float cImg = 2.0f * real[j] * img[j] + y[blockY*BLOCK_SIZE + i];
						float cReal = r2 - i2 + x[blockX*BLOCK_SIZE + j];

						isDone[j] = isDone[j] || outOfBound;
						doneCounter += isDone[j];

						data[index] = l * (!isDone[j]) + data[index] * isDone[j];

						img[j] = cImg * (!isDone[j]) + img[j] * isDone[j];
						real[j] = cReal * (!isDone[j]) + real[j] * isDone[j];
					}
				}
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
	memset(real, 0, BLOCK_SIZE * sizeof(float));
	memset(img, 0, BLOCK_SIZE * sizeof(float));
	memset(isDone, 0, BLOCK_SIZE * sizeof(bool));
	for(int i = 0; i < width; i++)
	{
		x[i] = x_start + i * dx;
	}
	for(int i = 0; i < height; i++)
	{
		y[i] = y_start + i * dy;
	}
}