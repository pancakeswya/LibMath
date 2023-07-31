#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#define S21_PI 3.14159265358979323846264338327950288
#define S21_EPS 1e-17
#define S21_EXP 2.7182818284590452353602874713526624
#define S21_E 2.71828182845904523536028747135266249
#define S21_DBL_MAX 1.7976931348623157e308L
#define S21_DBL_MIN 2.2250738585072014e-308L
#define S21_INF (1.0 / 0.0)
#define S21_NAN (0.0 / 0.0)

int s21_abs(int);
long double s21_fabs(double);
long double s21_sin(double);
long double s21_cos(double);
long double s21_tan(double);
long double s21_atan(double);
long double s21_acos(double);
long double s21_sqrt(double);
long double s21_asin(double);
long double s21_exp(double);
long double s21_log(double);
long double s21_fmod(double, double);
long double s21_pow(double, double);
long double s21_floor(double);
long double s21_ceil(double);

#endif  // SRC_S21_MATH_H_
