/**
 * @file BatchMandelCalculator.h
 * @author David Sladký (xsladk07@stud.fit.vutbr.cz)
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 2022-11-03
 */
#ifndef BATCHMANDELCALCULATOR_H
#define BATCHMANDELCALCULATOR_H

#include <BaseMandelCalculator.h>

class BatchMandelCalculator : public BaseMandelCalculator
{
public:
    BatchMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~BatchMandelCalculator();
    int * calculateMandelbrot();

private:
    int* data;
    float* real;
    float* img;
    float* y; 
    float* x; 
    float* oneBatchData;

    const size_t ALIGMENT = 64;
    const int BATCH_SIZE = 256;

    void initData();
};

#endif