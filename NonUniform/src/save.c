#include <errno.h> // errno, EEXIST
#include <stdbool.h> // false
#include <stdio.h> // snprintf
#include <string.h> // strlen
#include <sys/stat.h> // mode_t, S_IRWXU, S_IRWXG, S_IRWXO
#include "domain.h"
#include "flow_field.h"
#include "logger.h"
#include "memory.h"
#include "snpyio.h"
#include "./save.h"

#define ROOT_DIRECTORY "output/save/"

#define NDIMS 2

static int concat_dir_name(
    const size_t id,
    char ** const dir_name
) {
  const char prefix[] = ROOT_DIRECTORY;
  const int ndigits = 10;
  const int nchars = strlen(prefix) + ndigits + 1;
  *dir_name = memory_alloc(nchars, sizeof(char));
  (*dir_name)[nchars - 1] = '\0';
  if (nchars - 1 != snprintf(*dir_name, nchars, "%s%0*zu", prefix, ndigits, id)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

static int create_directory(
    const char dir_name[]
) {
  errno = 0;
  if (0 != mkdir(dir_name, S_IRWXU | S_IRWXG | S_IRWXO)) {
    perror(dir_name);
    // treat EEXIST as expected to override present files
    if (EEXIST != errno) {
      goto abort;
    }
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

static int write_npy_file(
    const char dir_name[],
    const char dset_name[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const size_t size,
    const void * data
){
  int error_code = 0;
  char * file_name = NULL;
  FILE * fp = NULL;
  size_t header_size = 0;
  // assign file_name
  {
    const char slash[] = {"/"};
    const char suffix[] = {".npy"};
    const int nchars =
      + strlen( dir_name)
      + strlen(    slash)
      + strlen(dset_name)
      + strlen(   suffix)
      + 1;
    file_name = memory_alloc(nchars, sizeof(char));
    file_name[nchars - 1] = '\0';
    if (nchars - 1 != snprintf(file_name, nchars, "%s%s%s%s", dir_name, slash, dset_name, suffix)) {
      error_code = 1;
      goto abort;
    }
  }
  // write npy header
  {
    errno = 0;
    fp = fopen(file_name, "w");
    if (NULL == fp) {
      perror(file_name);
      error_code = 1;
      goto abort;
    }
    size_t header_size = 0;
    if (0 != snpyio_w_header(ndims, shape, dtype, false, fp, &header_size)) {
      error_code = 1;
      goto abort;
    }
    fclose(fp);
  }
  // append main data
  {
    errno = 0;
    fp = fopen(file_name, "a");
    if (NULL == fp) {
      perror(file_name);
      error_code = 1;
      goto abort;
    }
    if (0 != fseek(fp, (long)header_size, SEEK_SET)) {
      error_code = 1;
      goto abort;
    }
    size_t nitems = 1;
    for (size_t dim = 0; dim < ndims; dim++) {
      nitems *= shape[dim];
    }
    if (nitems != fwrite(data, size, nitems, fp)) {
      error_code = 1;
      goto abort;
    }
  }
abort:
  memory_free(file_name);
  if (NULL != fp) {
    fclose(fp);
  }
  if (0 != error_code) {
    LOGGER_FAILURE();
  }
  return error_code;
}

int save(
    const size_t id,
    const size_t step,
    const double time,
    const domain_t * const domain,
    const flow_field_t * const flow_field
) {
  int error_code = 0;
  char * dir_name = NULL;
  if (0 != concat_dir_name(id, &dir_name)) {
    error_code = 1;
    goto abort;
  }
  if (0 != create_directory(dir_name)) {
    error_code = 1;
    goto abort;
  }
  write_npy_file(dir_name, "step", 0, NULL, "'<u8'", sizeof(size_t), &step);
  write_npy_file(dir_name, "time", 0, NULL, "'<f8'", sizeof(double), &time);
  write_npy_file(dir_name, "ux", NDIMS, (size_t [NDIMS]){NY + 2, NX + 2}, "'<f8'", sizeof(double), &flow_field->ux->buffer);
  write_npy_file(dir_name, "uy", NDIMS, (size_t [NDIMS]){NY + 2, NX + 2}, "'<f8'", sizeof(double), &flow_field->uy->buffer);
  write_npy_file(dir_name, "p", NDIMS, (size_t [NDIMS]){NY + 2, NX + 2}, "'<f8'", sizeof(double), &flow_field->p->buffer);
  write_npy_file(dir_name, "temperature", NDIMS, (size_t [NDIMS]){NY + 2, NX + 2}, "'<f8'", sizeof(double), &flow_field->temperature->buffer);
  write_npy_file(dir_name, "salinity", NDIMS, (size_t [NDIMS]){NY + 2, NX + 2}, "'<f8'", sizeof(double), &flow_field->salinity->buffer);
  if (0 == id) {
    write_npy_file(ROOT_DIRECTORY, "x_centers", 1, (size_t [NDIMS]){NX + 2}, "'<f8'", sizeof(double), domain->x_coordinate.centers);
    write_npy_file(ROOT_DIRECTORY, "x_faces", 1, (size_t [NDIMS]){NX + 2}, "'<f8'", sizeof(double), domain->x_coordinate.faces);
    write_npy_file(ROOT_DIRECTORY, "y_centers", 1, (size_t [NDIMS]){NY + 2}, "'<f8'", sizeof(double), domain->y_coordinate.centers);
    write_npy_file(ROOT_DIRECTORY, "y_faces", 1, (size_t [NDIMS]){NY + 2}, "'<f8'", sizeof(double), domain->y_coordinate.faces);
  }
abort:
  memory_free(dir_name);
  if (0 != error_code) {
    LOGGER_FAILURE();
  }
  return error_code;
}

