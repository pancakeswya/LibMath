#include "s21_math.h"

int s21_abs(int x) { return x < 0 ? -x : x; }

long double s21_fabs(double x) { return x < 0 ? -x : x; }

long double s21_sin(double x) {
  double result, elem;
  if (s21_fabs(x) == S21_INF) {
    result = -S21_NAN;
  } else {
    x = (x > 2 * S21_PI || x < -2 * S21_PI) ? s21_fmod(x, 2 * S21_PI) : x;
    result = elem = x;
    for (int i = 1; s21_fabs(elem) > S21_EPS; i++) {
      elem *= x * x / ((2.0 * i) * ((2.0 * i) + 1));
      result += i % 2 ? -elem : elem;
    }
  }
  return result;
}

long double s21_cos(double x) {
  return s21_fabs(x) == S21_INF ? -S21_NAN : s21_sin(S21_PI / 2. - x);
}

long double s21_tan(double x) { return s21_sin(x) / s21_cos(x); }

long double s21_atan(double x) {
  long double res, ch, elem;
  int flag_span = 0, sign = (x < 0) ? 0 : 1;
  x = s21_fabs(x);
  res = ch = elem = x;
  if (x > 1) {
    flag_span = 1;
    ch = res = x = 1 / x;
  }
  if (x == 1) {
    res = S21_PI / 4;
  } else {
    for (int n = 1; s21_fabs(elem) > S21_EPS; n++) {
      ch *= x * x;
      elem = ch / (n * 2 + 1);
      res += n % 2 ? -elem : elem;
    }
  }
  if (flag_span) {
    res = S21_PI / 2 - res;
  }
  if (!sign) {
    res = -res;
  }
  return res;
}

long double s21_sqrt(double x) {
  long double rez = 4, y = 0;
  if (x == S21_INF) {
    rez = S21_INF;
  } else if (x == -S21_INF) {
    rez = S21_NAN;
  } else {
    while (s21_fabs(rez - y) > S21_EPS) {
      if (x < 0) {
        rez = S21_NAN;
        break;
      }
      y = rez;
      rez = (y + x / y) / 2;
    }
  }
  if (x == 0) {
    rez = 0;
  }
  return rez;
}

long double s21_asin(double x) {
  long double result;
  int sign = (x < 0) ? 0 : 1;
  x = s21_fabs(x);
  if (x == 1) {
    result = S21_PI / 2;
  } else if (x < 1 && x >= 0) {
    result = s21_atan(x / s21_sqrt(1 - x * x));
  } else {
    result = S21_NAN;
    x == S21_INF ? sign = 0 : sign;
  }
  if (!sign) {
    result = - result;
  }
  return result;
}

long double s21_acos(double x) {
  long double res;
  if ((x <= 1) && (x >= -1)) {
    res = S21_PI / 2 - s21_asin(x);
  } else if (s21_fabs(x) == S21_INF) {
    res = -(S21_NAN);
  } else {
    res = S21_NAN;
  }
  return res;
}

long double s21_exp(double x) {
  long double res = 1, y = 1, n = 1;
  int neg = (x < 0);
  x = s21_fabs(x);
  while (s21_fabs((double)res) > S21_EPS) {
    if (y > S21_DBL_MAX) {
      y = S21_INF;
      break;
    }
    res *= x / n++;
    y += res;
  }
  y = (neg) ? y > S21_DBL_MAX ? 0 : 1. / y : y;
  return y;
}

long double s21_log(double x) {
  long double res;
  if (x == S21_INF) {
    res = S21_INF;
  } else if (x == -S21_INF) {
    res = S21_NAN;
  } else {
    unsigned a = 0, d = 1;
    long double b = 0, c = (x < 1) ? 1 / x : x, e = b, f;
    if (x > 0) {
      for (; (c /= S21_EXP) > 1; ++a)
        ;
      c = 2 * (1 / (c * S21_EXP - 1)) + 1;
      f = c * c;
      for (c /= 2; b += 1 / (d * c), b - e; d += 2) {
        c *= f;
        e = b;
      }
    } else
      b = (x == 0) / 0.;
    res = a + b;
  }
  if ((s21_fabs(x) != S21_INF) && (x < 1)) {
    res = -res;
  }
  return res;
}

long double s21_fmod(double x, double y) {
  long double rem;
  if (x == S21_INF || y == 0) {
    rem = S21_NAN;
  } else if (y == S21_INF || y == -S21_INF) {
    rem = (long long)x;
  } else {
    long double dx = x;
    long double dy = y;
    long long int div = dx / dy;
    dx -= div * dy;
    rem = dx;
  }
  return rem;
}

long double s21_pow(double base, double exp) {
  long double power;
  if ((long int)exp == exp) {
    if (exp > 0) {
      power = base;
      for (long int i = 0; i < (long int)exp - 1; i++) {
        power *= base;
      }
    } else if (exp == 0) {
      power = 1;
    } else {
      power = 1 / base;
      for (long int i = 0; i < (long int)exp * (-1) - 1; i++) {
        power /= base;
      }
    }
  } else if (base < 0) {
    if (s21_fabs(exp) == S21_INF) {
      if (s21_fabs(base) < 1) {
        power = 0;
      } else if (s21_fabs(base) == 1) {
        power = 1;
      } else {
        power = (exp == -S21_INF) ? 0 : S21_INF;
      }
    } else {
      power = -S21_NAN;
    }
  } else if (base == 0) {
    power = (exp == 0) ? 1 : 0;
  } else if (base == 1) {
    power = 1;
  } else {
    power = s21_exp(exp * (double)s21_log(base));
  }
  return power;
}

long double s21_floor(double x) {
  long double res = x;
  if (s21_fabs(x) == S21_INF) {
    res = x;
  } else if (x == S21_NAN) {
    res = S21_NAN;
  } else {
    if (res >= 0) {
      res -= s21_fmod(x, 1);
    } else if (s21_fabs(s21_fmod(x, 1)) > S21_EPS) {
      res -= 1;
      res -= s21_fmod(x, 1);
    }
  }
  return res;
}

long double s21_ceil(double x) {
  long double res = x;
  if (s21_fabs(x) == S21_INF) {
    res = x;
  } else if (x == S21_NAN) {
    res = S21_NAN;
  } else if (s21_fabs(s21_fmod(x, 1)) > 0) {
    if (x > 0.0) {
      res = (long int)(x + 1.0);
    } else {
      res = (long int)x;
    }
  }
  return res;
}
