
#ifndef CONTINUOUS_H
#define CONTINUOUS_H


#include <math.h>


double abs_control(int *sources, int *targets,
        const int num_links, double *expression);
double difference_control(int *sources, int *targets,
        const int num_links, double *expression);
double abs_difference_control(int *sources, int *targets,
        const int num_links, double *expression);
double functional_control(int *sources, int *targets, int *function,
        const int num_links, double *expression);
double functional_comparison(int *sources, int *targets, int *function,
        const int num_links, double *rate);
double delayed_functional_control(int *sources, int *targets, int *function,
        const int num_links, double *expression, double *delta_expression);


#endif // CONTINUOUS_H
