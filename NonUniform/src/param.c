#include <math.h>
#include "param.h"

int setup_parameters(
    const non_dimensional_numbers_t * const non_dimensional_numbers,
    parameters_t * const parameters
) {
  parameters->density_ratio = non_dimensional_numbers->density_ratio;
  parameters->fluid_diffusivity = sqrt(
      non_dimensional_numbers->prandtl / non_dimensional_numbers->rayleigh
  );
  parameters->salinity_diffusivity = 1. / non_dimensional_numbers->lewis * sqrt(
      1. / non_dimensional_numbers->prandtl / non_dimensional_numbers->rayleigh
  );
  parameters->temperature_diffusivity = sqrt(
      1. / non_dimensional_numbers->prandtl / non_dimensional_numbers->rayleigh
  );
  return 0;
}

