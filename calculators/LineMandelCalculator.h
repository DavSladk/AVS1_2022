/**
 * @file LineMandelCalculator.h
 * @author David Sladký (xsladk07@stud.fit.vutbr.cz)
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 2022-11-03
 */

#include <BaseMandelCalculator.h>

class LineMandelCalculator : public BaseMandelCalculator
{
public:
    LineMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~LineMandelCalculator();
    int *calculateMandelbrot();

private:
    int* data;
    float* real;
    float* img;
    float* y; 
    float* x; 
    float* oneLineData;

    const size_t ALIGMENT = 64;

    void initData();
};