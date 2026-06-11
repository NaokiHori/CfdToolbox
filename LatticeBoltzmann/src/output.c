#include <math.h>
#include <stdlib.h>
#include "./output.h"
#include "./output/snpyio.h"

int output_distribution_function(
    const char * file_name,
    const df_t * const df
) {
  double * const buf = malloc(NDIRS * NX * (NY + 2) * sizeof(double));
  for (size_t cnt = 0, j = 0; j <= NY + 1; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double * const dfij = (*df)[j][i];
      for (size_t n = 0; n < NDIRS; n++) {
        buf[cnt++] = dfij[n];
      }
    }
  }
  FILE * const fp = fopen(file_name, "w");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  size_t header_size = 0;
  const size_t ndims[NDIMS + 1] = {NY + 2, NX, NDIRS};
  if (0 != snpyio_w_header(NDIMS + 1, ndims, "'<f8'", false, fp, &header_size)) {
    return 1;
  }
  const size_t nitems = NDIRS * NX * (NY + 2);
  if (nitems != fwrite(buf, sizeof(double), nitems, fp)) {
    return 1;
  }
  free(buf);
  fclose(fp);
  return 0;
}

