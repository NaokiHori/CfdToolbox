#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include "memory.h"
#include "./output.h"
#include "./output/snpyio.h"

#define FILE_NAME_MAX_SIZE 256
#define ROOT_DIRECTORY "output/"
#define NDIMS 2

static int output_vector(
    const char vector_name[],
    const size_t nitems,
    const double * const vector
) {
  char file_name[FILE_NAME_MAX_SIZE] = {'\0'};
  snprintf(file_name, FILE_NAME_MAX_SIZE - 1, ROOT_DIRECTORY "%s.npy", vector_name);
  errno = 0;
  FILE * const fp = fopen(file_name, "w");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  size_t header_size = 0;
  if (0 != snpyio_w_header(1, &nitems, "'<f8'", false, fp, &header_size)) {
    fclose(fp);
    return 1;
  }
  fwrite(vector, nitems, sizeof(double), fp);
  fclose(fp);
  return 0;
}

static int output_array(
    const char array_name[],
    const size_t counter,
    const size_t nitems[NDIMS],
    double ** const array
) {
  char file_name[FILE_NAME_MAX_SIZE] = {'\0'};
  snprintf(file_name, FILE_NAME_MAX_SIZE - 1, ROOT_DIRECTORY "%s_%010zu.npy", array_name, counter);
  errno = 0;
  FILE * const fp = fopen(file_name, "w");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  size_t header_size = 0;
  if (0 != snpyio_w_header(NDIMS, nitems, "'<f8'", false, fp, &header_size)) {
    fclose(fp);
    return 1;
  }
  fwrite(&array[0][0], nitems[0] * nitems[1], sizeof(double), fp);
  fclose(fp);
  return 0;
}

int output(
    const domain_t * const domain,
    double ** const concentration
) {
  static size_t counter = 0;
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dy = domain->dy;
  if (0 == counter) {
    double * const yc = memory_alloc(ny, sizeof(double));
    for (size_t j = 0; j < ny; j++) {
      yc[j] = 0.5 * (2 * j + 1) * dy;
    }
    output_vector("xc", nx + 2, domain->xc);
    output_vector("yc", ny, yc);
    memory_free(yc);
  }
  output_array("concentration", counter, (size_t [NDIMS]){nx + 2, ny + 2}, concentration);
  counter += 1;
  return 0;
}

