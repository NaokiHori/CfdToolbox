#include <errno.h>
#include <stdio.h>
#include "./output.h"
#include "./output/snpyio.h"

static int compute_total_nitems(
    const size_t ndims,
    const size_t * const nitems,
    size_t * const total_nitems
) {
  *total_nitems = 1;
  for (size_t dim = 0; dim < ndims; dim++) {
    *total_nitems *= nitems[dim];
  }
  return 0;
}

int output(
    const char file_name[],
    const size_t ndims,
    const size_t * const nitems,
    const double * const buffer
) {
  errno = 0;
  FILE * const fp = fopen(file_name, "w");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  size_t header_size = 0;
  if (0 != snpyio_w_header(ndims, nitems, "'<f8'", false, fp, &header_size)) {
    return 1;
  }
  size_t total_nitems = 0;
  compute_total_nitems(ndims, nitems, &total_nitems);
  fwrite(buffer, sizeof(double), total_nitems, fp);
  fclose(fp);
  return 0;
}
