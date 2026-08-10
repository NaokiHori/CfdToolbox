#include "domain.h"
#include "load_grid.h"
#include "memory.h"
#include "param.h"

// grid-streching indicator
// 0: extremely stretched
// 1: uniform
static int build_coordinate(
    const double length,
    const size_t nitems,
    const char centers_file[],
    const char faces_file[],
    const char center_scale_factors_file[],
    const char face_scale_factors_file[],
    coordinate_t * const coordinate
) {
  coordinate->length = length;
  // coordinates
  double ** const centers = &coordinate->centers;
  double ** const faces = &coordinate->faces;
  *centers = memory_alloc(nitems + 2, sizeof(double));
  *faces = memory_alloc(nitems + 2, sizeof(double));
  if (0 != try_load_npy(
      centers_file,
      1,
      (size_t [1]){nitems + 2},
      "'<f8'",
      false,
      sizeof(double),
      *centers
  )) {
    return 1;
  }
  if (0 != try_load_npy(
      faces_file,
      1,
      (size_t [1]){nitems + 2},
      "'<f8'",
      false,
      sizeof(double),
      *faces
  )) {
    return 1;
  }
  // scale factors
  double ** const center_scale_factors = &coordinate->center_scale_factors;
  double ** const face_scale_factors = &coordinate->face_scale_factors;
  *center_scale_factors = memory_alloc(nitems + 2, sizeof(double));
  *face_scale_factors = memory_alloc(nitems + 2, sizeof(double));
  if (0 != try_load_npy(
      center_scale_factors_file,
      1,
      (size_t [1]){nitems + 2},
      "'<f8'",
      false,
      sizeof(double),
      *center_scale_factors
  )) {
    return 1;
  }
  if (0 != try_load_npy(
      face_scale_factors_file,
      1,
      (size_t [1]){nitems + 2},
      "'<f8'",
      false,
      sizeof(double),
      *face_scale_factors
  )) {
    return 1;
  }
  return 0;
}

static int coordinate_finalize(
    coordinate_t * const coordinate
) {
  memory_free(coordinate->faces);
  memory_free(coordinate->centers);
  memory_free(coordinate->face_scale_factors);
  memory_free(coordinate->center_scale_factors);
  return 0;
}

int domain_init(
    domain_t * const domain
) {
  double lengths[NDIMS] = {0., 0.};
  size_t nitems[NDIMS] = {0, 0};
  if (0 != try_load_npy(
      GRID_DIRECTORY "lengths.npy",
      1,
      (size_t [1]){NDIMS},
      "'<f8'",
      false,
      sizeof(double),
      lengths
  )) {
    return 1;
  }
  if (0 != try_load_npy(
      GRID_DIRECTORY "nitems.npy",
      1,
      (size_t [1]){NDIMS},
      "'<u8'",
      false,
      sizeof(size_t),
      nitems
  )) {
    return 1;
  }
  if (0 != build_coordinate(
      lengths[0],
      nitems[0],
      GRID_DIRECTORY "x_centers.npy",
      GRID_DIRECTORY "x_faces.npy",
      GRID_DIRECTORY "x_center_scale_factors.npy",
      GRID_DIRECTORY "x_face_scale_factors.npy",
      &domain->x_coordinate
  )) {
    return 1;
  }
  if (0 != build_coordinate(
      lengths[1],
      nitems[1],
      GRID_DIRECTORY "y_centers.npy",
      GRID_DIRECTORY "y_faces.npy",
      GRID_DIRECTORY "y_center_scale_factors.npy",
      GRID_DIRECTORY "y_face_scale_factors.npy",
      &domain->y_coordinate
  )) {
    return 1;
  }
  domain->nx = nitems[0];
  domain->ny = nitems[1];
  return 0;
}

int domain_finalize(
    domain_t * const domain
) {
  coordinate_finalize(&domain->x_coordinate);
  coordinate_finalize(&domain->y_coordinate);
  return 0;
}

