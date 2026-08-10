#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "load_grid.h"
#include "snpyio.h"

int try_load_npy(
    const char file_name[],
    const size_t ndims,
    const size_t * const shape,
    const char dtype[],
    const bool is_fortran_order,
    const size_t element_size,
    void * const buffer
) {
  errno = 0;
  FILE * const fp = fopen(file_name, "r");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  size_t obtained_ndims = 0;
  size_t * obtained_dims = NULL;
  char * obtained_dtype = NULL;
  bool obtained_is_fortran_order = false;
  size_t header_size = 0;
  if (0 != snpyio_r_header(
      &obtained_ndims,
      &obtained_dims,
      &obtained_dtype,
      &obtained_is_fortran_order,
      fp,
      &header_size
  )) {
    printf("file %s: failed to load numpy header\n", file_name);
    goto abort;
  }
  // assert header
  if (ndims != obtained_ndims) {
    printf("file %s: unexpected ndims, expected %zu obtained %zu\n", file_name, ndims, obtained_ndims);
    goto abort;
  }
  for (size_t n = 0; n < ndims; n++) {
    if (shape[n] != obtained_dims[n]) {
      printf(
          "file %s: unexpected shape[%zu], expected %zu obtained %zu\n",
          file_name,
          n,
          shape[n],
          obtained_dims[n]
      );
      goto abort;
    }
  }
  if (0 != strcmp(dtype, obtained_dtype)) {
    printf("file %s: unexpected dtype, expected %s obtained %s\n", file_name, dtype, obtained_dtype);
    goto abort;
  }
  if (is_fortran_order != obtained_is_fortran_order) {
    printf(
        "file %s: unexpected is_fortran_order, expected %s obtained %s\n",
        file_name,
        is_fortran_order ? "true" : "false",
        obtained_is_fortran_order ? "true" : "false"
    );
    goto abort;
  }
  // load buffer
  size_t nitems = 1;
  for (size_t n = 0; n < ndims; n++) {
    nitems *= obtained_dims[n];
  }
  fread(buffer, element_size, nitems, fp);
  fclose(fp);
  free(obtained_dims);
  free(obtained_dtype);
  return 0;
abort:
  fclose(fp);
  free(obtained_dims);
  free(obtained_dtype);
  return 1;
}

