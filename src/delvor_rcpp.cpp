#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include "deltri3.h"

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

static inline size_t next_halfedge(size_t e) {
  return (e % 3 == 2) ? e - 2 : e + 1;
}

static inline size_t prev_halfedge(size_t e) {
  return (e % 3 == 0) ? e + 2 : e - 1;
}

static inline size_t triangle_of_edge(size_t e) {
  return e / 3;
}

struct Circumcenter {
  double x;
  double y;
};

static inline Circumcenter circumcenter_pts(
    double ax, double ay,
    double bx, double by,
    double cx, double cy) {
  
  const double dx = bx - ax;
  const double dy = by - ay;
  const double ex = cx - ax;
  const double ey = cy - ay;
  const double bl = dx * dx + dy * dy;
  const double cl = ex * ex + ey * ey;
  const double d  = dx * ey - dy * ex;
  
  Circumcenter out;
  out.x = ax + (ey * bl - dy * cl) * 0.5 / d;
  out.y = ay + (dx * cl - ex * bl) * 0.5 / d;
  return out;
}

// point-in-polygon test for convex hull sign choice
static inline bool point_in_polygon(
    double px, double py,
    const std::vector<double>& hx,
    const std::vector<double>& hy) {
  
  bool inside = false;
  const size_t n = hx.size();
  if (n < 3) return false;
  
  for (size_t i = 0, j = n - 1; i < n; j = i++) {
    const bool intersect =
      ((hy[i] > py) != (hy[j] > py)) &&
      (px < (hx[j] - hx[i]) * (py - hy[i]) / (hy[j] - hy[i] + 1e-300) + hx[i]);
    
    if (intersect) inside = !inside;
  }
  
  return inside;
}

// closer port of alphahull::dummycoor()
static inline Circumcenter dummy_point_exact(
    double x1, double y1, double x2, double y2,
    double cx, double cy,
    double away,
    const std::vector<double>& hull_x,
    const std::vector<double>& hull_y) {
  
  const double eps = 1e-5;
  
  // v <- l2 - l1 ; then c(v[2], -v[1])
  double vx = x2 - x1;
  double vy = y2 - y1;
  double nx = vy;
  double ny = -vx;
  
  // alphahull::dummycoor divides by sum(v^2), not sqrt(sum(v^2))
  const double norm = nx * nx + ny * ny;
  if (norm > 0.0) {
    nx /= norm;
    ny /= norm;
  }
  
  const double mx = 0.5 * (x1 + x2);
  const double my = 0.5 * (y1 + y2);
  
  const double tx = mx + eps * nx;
  const double ty = my + eps * ny;
  
  Circumcenter out;
  if (point_in_polygon(tx, ty, hull_x, hull_y)) {
    out.x = cx - away * nx;
    out.y = cy - away * ny;
  } else {
    out.x = cx + away * nx;
    out.y = cy + away * ny;
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::List delvor_rcpp(Rcpp::NumericMatrix x) {
  const size_t n = x.nrow();
  if (x.ncol() != 2) stop("x must be an n x 2 matrix");
  if (n < 3) stop("At least three points are required");
  
  std::vector<double> coords(2 * n);
  double xmin = x(0, 0), xmax = x(0, 0), ymin = x(0, 1), ymax = x(0, 1);
  
  for (size_t i = 0; i < n; ++i) {
    coords[2 * i]     = x(i, 0);
    coords[2 * i + 1] = x(i, 1);
    xmin = std::min(xmin, x(i, 0));
    xmax = std::max(xmax, x(i, 0));
    ymin = std::min(ymin, x(i, 1));
    ymax = std::max(ymax, x(i, 1));
  }
  
  deltri::deltri_cpp del(coords);
  const auto& tri  = del.triangles;
  const auto& half = del.halfedges;
  
  const size_t nt = tri.size() / 3;
  std::vector<Circumcenter> cc(nt);
  
  for (size_t t = 0; t < nt; ++t) {
    const size_t i0 = tri[3 * t];
    const size_t i1 = tri[3 * t + 1];
    const size_t i2 = tri[3 * t + 2];
    
    cc[t] = circumcenter_pts(
      coords[2 * i0], coords[2 * i0 + 1],
                            coords[2 * i1], coords[2 * i1 + 1],
                                                  coords[2 * i2], coords[2 * i2 + 1]
    );
  }
  
  // Rebuild ordered convex hull using hull_start / hull_next
  std::vector<size_t> hull_ids;
  {
    size_t start = del.hull_start;
    if (start != deltri::INVALID_INDEX) {
      size_t v = start;
      do {
        hull_ids.push_back(v);
        v = del.hull_next[v];
      } while (v != start && v != deltri::INVALID_INDEX);
    }
  }
  
  std::vector<double> hull_x;
  std::vector<double> hull_y;
  hull_x.reserve(hull_ids.size() + 1);
  hull_y.reserve(hull_ids.size() + 1);
  
  for (size_t k = 0; k < hull_ids.size(); ++k) {
    const size_t id = hull_ids[k];
    hull_x.push_back(coords[2 * id]);
    hull_y.push_back(coords[2 * id + 1]);
  }
  
  // close polygon
  if (!hull_x.empty()) {
    hull_x.push_back(hull_x.front());
    hull_y.push_back(hull_y.front());
  }
  
  // Hull successor map
  std::unordered_map<size_t, size_t> hull_succ;
  for (size_t k = 0; k < hull_ids.size(); ++k) {
    hull_succ[hull_ids[k]] = hull_ids[(k + 1) % hull_ids.size()];
  }
  
  const double away = std::max(xmax - xmin, ymax - ymin);
  
  struct Row {
    size_t ind1, ind2;
    double x1, y1, x2, y2;
    double mx1, my1, mx2, my2;
    int bp1, bp2;
  };
  
  std::vector<Row> rows;
  rows.reserve(tri.size() / 2);
  
  for (size_t e = 0; e < tri.size(); ++e) {
    const size_t opp = half[e];
    
    // Keep only one representative for each undirected edge
    if (opp != deltri::INVALID_INDEX && e > opp) continue;
    
    size_t a = tri[e];
    size_t b = tri[next_halfedge(e)];
    
    // Force hull-edge orientation to follow the ordered hull cycle
    if (opp == deltri::INVALID_INDEX) {
      auto it = hull_succ.find(a);
      if (it != hull_succ.end() && it->second == b) {
        // already correctly oriented
      } else {
        auto it2 = hull_succ.find(b);
        if (it2 != hull_succ.end() && it2->second == a) {
          std::swap(a, b);
        }
      }
    }
    
    const size_t t1 = triangle_of_edge(e);
    const Circumcenter c1 = cc[t1];
    
    int bp1 = 0;
    int bp2 = 0;
    Circumcenter c2;
    
    if (opp == deltri::INVALID_INDEX) {
      bp2 = 1;
      c2 = dummy_point_exact(
        coords[2 * a], coords[2 * a + 1],
                             coords[2 * b], coords[2 * b + 1],
                                                  c1.x, c1.y,
                                                  away,
                                                  hull_x, hull_y
      );
    } else {
      const size_t t2 = triangle_of_edge(opp);
      c2 = cc[t2];
    }
    
    Row r;
    r.ind1 = a + 1; // R is 1-based
    r.ind2 = b + 1;
    r.x1   = coords[2 * a];
    r.y1   = coords[2 * a + 1];
    r.x2   = coords[2 * b];
    r.y2   = coords[2 * b + 1];
    r.mx1  = c1.x;
    r.my1  = c1.y;
    r.mx2  = c2.x;
    r.my2  = c2.y;
    r.bp1  = bp1;
    r.bp2  = bp2;
    
    rows.push_back(r);
  }
  
  NumericMatrix mesh(rows.size(), 12);
  colnames(mesh) = CharacterVector::create(
    "ind1", "ind2", "x1", "y1", "x2", "y2",
    "mx1", "my1", "mx2", "my2", "bp1", "bp2"
  );
  
  for (size_t i = 0; i < rows.size(); ++i) {
    mesh(i, 0)  = rows[i].ind1;
    mesh(i, 1)  = rows[i].ind2;
    mesh(i, 2)  = rows[i].x1;
    mesh(i, 3)  = rows[i].y1;
    mesh(i, 4)  = rows[i].x2;
    mesh(i, 5)  = rows[i].y2;
    mesh(i, 6)  = rows[i].mx1;
    mesh(i, 7)  = rows[i].my1;
    mesh(i, 8)  = rows[i].mx2;
    mesh(i, 9)  = rows[i].my2;
    mesh(i,10)  = rows[i].bp1;
    mesh(i,11)  = rows[i].bp2;
  }
  
  return List::create(
    _["mesh"] = mesh,
    _["x"] = x,
    _["tri.obj"] = R_NilValue
  );
}