
#include "continuous.h"


double
abs_control(int *sources, int *targets, const int num_links, double *expression)
{
    double sum = 0.0;
    int i = 0;
    for (i = 0; i < num_links; ++i) {
        sum += 1.0 - fabs(expression[sources[i]] - expression[targets[i]]);
    }
    return sum / (double)num_links;
}

double
difference_control(int *sources, int *targets, const int num_links,
        double *expression)
{
    double sum = 0.0;
    int i = 0;
    for (i = 0; i < num_links; ++i) {
        sum += expression[sources[i]] - expression[targets[i]];
    }
    return sum / (double)num_links;
}

double
abs_difference_control(int *sources, int *targets, const int num_links,
        double *expression)
{
    double sum = 0.0;
    int i = 0;
    for (i = 0; i < num_links; ++i) {
        sum += fabs(expression[sources[i]] - expression[targets[i]]);
    }
    return sum / (double)num_links;
}

double
functional_control(int *sources, int *targets, int *function, const int num_links,
        double *expression)
{
    double sum = 0.0;
    int i = 0;
    for (i = 0; i < num_links; ++i) {
        if (function[i] == 1) {
            sum += 1.0 - fabs(expression[sources[i]] - expression[targets[i]]);
        }
        else if (function[i] == -1) {
            sum += fabs(expression[sources[i]] - expression[targets[i]]);
        }
    }
    return sum / (double)num_links;
}

double
functional_comparison(int *sources, int *targets, int *function, const int num_links,
        double *rate)
{
    double sum = 0.0;
    int i = 0;
    for (i = 0; i < num_links; ++i) {
        if (function[i] == 1) {
            sum += (rate[sources[i]] > 0) ? (double)(rate[targets[i]] >= 0.0) : (double)(rate[targets[i]] <= 0.0);
        }
        else if (function[i] == -1) {
            sum += (rate[sources[i]] > 0) ? (double)(rate[targets[i]] <= 0.0) : (double)(rate[targets[i]] >= 0.0);
        }
    }
    return sum / (double)num_links;
}

double
delayed_functional_control(int *sources, int *targets, int *function, const int num_links,
        double *expression, double *delta_expression)
{
    double sum = 0.0;
    int i = 0;
    for (i = 0; i < num_links; ++i) {
        if (function[i] == 1) {
            sum += 1.0 - fabs(expression[sources[i]] - delta_expression[targets[i]]);
        }
        else if (function[i] == -1) {
            sum += fabs(expression[sources[i]] - delta_expression[targets[i]]);
        }
    }
    return sum / (double)num_links;
}

