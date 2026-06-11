#include "time_marcher.h"

// three-step scheme by Williamson, JCP, 1980
const runge_kutta_coefs_t RUNGE_KUTTA_COEFS = {
  .alpha = {+ 1920. / 5760., + 2400. / 5760., + 1440. / 5760.},
  .beta  = {     0. / 5760., - 3200. / 5760., - 6885. / 5760.},
  .gamma = {+ 1920. / 5760., + 5400. / 5760., + 3072. / 5760.},
};

