/**
 * @file LineMandelCalculator.h
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */

#include <BaseMandelCalculator.h>

#define REDUCE

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
    
    #ifdef REDUCE
    bool* done;
    #endif

    const size_t ALIGMENT = 64;

    void initData();
};