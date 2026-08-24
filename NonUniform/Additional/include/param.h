#if !defined(PARAM_H)
#define PARAM_H

#define PI 3.141592653598793238

typedef struct {
  const double density_ratio;
  const double lewis;
  const double prandtl;
  const double rayleigh;
} non_dimensional_numbers_t;

typedef struct {
  double density_ratio;
  double fluid_diffusivity;
  double salinity_diffusivity;
  double temperature_diffusivity;
} parameters_t;

extern int setup_parameters(
    const non_dimensional_numbers_t * const non_dimensional_numbers,
    parameters_t * const parameters
);

#endif // PARAM_H
