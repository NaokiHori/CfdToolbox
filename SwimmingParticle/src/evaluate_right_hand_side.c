#include "./evaluate_right_hand_side.h"

int evaluate_right_hand_side(
    const double peclet,
    const domain_t * const domain,
    const size_t runge_kutta_step,
    double ** const concentration,
    double ** const stream_function,
    runge_kutta_buffer_t * const runge_kutta_buffer
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const hxxf = domain->hxxf;
  const double * const hxxc = domain->hxxc;
  const double * const hyxf = domain->hyxf;
  const double * const hyxc = domain->hyxc;
  const double * const jdxf = domain->jdxf;
  const double * const jdxc = domain->jdxc;
  double ** const non_linear = runge_kutta_buffer->non_linear;
  double ** const linear = runge_kutta_buffer->increment;
  if (0 == runge_kutta_step) {
#pragma omp parallel for
    for (size_t i = 1; i <= nx; i++) {
      for (size_t j = 1; j <= ny; j++) {
        non_linear[i][j] = 0.;
      }
    }
  }
  const double rk_b = RUNGE_KUTTA_COEFS.beta[runge_kutta_step];
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t j = 1; j <= ny; j++) {
      // compute velocity by differentiating the stream function
      //   u1 = + 1 / h2 * dsdx2
      //   u2 = - 1 / h1 * dsdx1
      const double ux_xm = + 1. / hyxf[i    ] * (
          - stream_function[i    ][j    ]
          + stream_function[i    ][j + 1]
      );
      const double ux_xp = + 1. / hyxf[i + 1] * (
          - stream_function[i + 1][j    ]
          + stream_function[i + 1][j + 1]
      );
      const double uy_ym = - 1. / hxxc[i    ] * (
          - stream_function[i    ][j    ]
          + stream_function[i + 1][j    ]
      );
      const double uy_yp = - 1. / hxxc[i    ] * (
          - stream_function[i    ][j + 1]
          + stream_function[i + 1][j + 1]
      );
      // concentration gradient on general curvilinear coordinate system
      const double dc_xm = (
          - concentration[i - 1][j    ]
          + concentration[i    ][j    ]
      );
      const double dc_xp = (
          - concentration[i    ][j    ]
          + concentration[i + 1][j    ]
      );
      const double dc_ym = (
          - concentration[i    ][j - 1]
          + concentration[i    ][j    ]
      );
      const double dc_yp = (
          - concentration[i    ][j    ]
          + concentration[i    ][j + 1]
      );
      // advective contribution
      const double advx = - 1. / jdxc[i    ] * (
          + 0.5 * jdxf[i    ] / hxxf[i    ] * ux_xm * dc_xm
          + 0.5 * jdxf[i + 1] / hxxf[i + 1] * ux_xp * dc_xp
      );
      const double advy = - 1. / jdxc[i    ] * (
          + 0.5 * jdxc[i    ] / hyxc[i    ] * uy_ym * dc_ym
          + 0.5 * jdxc[i    ] / hyxc[i    ] * uy_yp * dc_yp
      );
      // diffusive contribution
      const double difx = 1. / peclet / jdxc[i    ] * (
          - jdxf[i    ] / hxxf[i    ] / hxxf[i    ] * dc_xm
          + jdxf[i + 1] / hxxf[i + 1] / hxxf[i + 1] * dc_xp
      );
      const double dify = 1. / peclet / jdxc[i    ] * (
          - jdxc[i    ] / hyxc[i    ] / hyxc[i    ] * dc_ym
          + jdxc[i    ] / hyxc[i    ] / hyxc[i    ] * dc_yp
      );
      // non-linear terms are appended to the previous information
      non_linear[i][j] = rk_b * non_linear[i][j] + advx + advy;
      // linear terms do not need past information
      linear[i][j] = difx + dify;
    }
  }
  return 0;
}

