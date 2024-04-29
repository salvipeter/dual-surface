#include <cmath>
#include <fstream>

#include <c0coons.hh>

using namespace Geometry;

struct Plane {
  Point3D p;  // any point on the plane
  Vector3D n; // unit normal
};

struct Line {
  Point3D p;  // any point on the line
  Vector3D d; // unit tangent
};

using SPoint = Point3D; // spherical point
// represents a plane with (r, theta, phi) s.t.:
//         r       is the distance to the origin
//   (theta, phi) are the spherical coordinates of the unit normal

SPoint planeToSPoint(const Plane &p) {
  auto r = p.p * p.n;
  auto theta = std::acos(p.n[2]), phi = 0.0;
  auto tmp = p.n;
  tmp[2] = 0;
  auto denom = tmp.norm();
  if (denom > 0)
    phi = std::acos(p.n[0] / tmp.norm()) * (p.n[1] < 0 ? -1 : 1);
  return { r, theta, phi };
}

Plane spointToPlane(const SPoint &s) {
  Plane p;
  p.n[0] = std::sin(s[1]) * std::cos(s[2]);
  p.n[1] = std::sin(s[1]) * std::sin(s[2]);
  p.n[2] = std::cos(s[1]);
  p.p = p.n * s[0];
  return p;
}

Line intersectPlanes(const Plane &p1, const Plane &p2) {
  Line l;
  l.p = p2.p - p1.n * ((p2.p - p1.p) * p1.n);
  l.d = p1.n ^ p2.n;
  auto denom = l.d.norm();
  if (denom > 0)
    l.d /= denom;
  return l;
}

// Midpoint between closest points when the lines do not intersect
Point3D intersectLines(const Line &l1, const Line &l2) {
  auto ap = l1.p, bp = l2.p;
  auto ad = l1.d, bd = l2.d;
  double a = ad * ad, b = ad * bd, c = bd * bd;
  double d = ad * (ap - bp), e = bd * (ap - bp);
  if (a * c - b * b < 1.0e-7)
    return ap;
  double D = a * c - b * b;
  double s = (b * e - c * d) / D;
  double t = (a * e - b * d) / D;
  return ((ap + ad * s) + (bp + bd * t)) / 2;
}


// Ribbon definition & import from file

class NormalFence : public CurveType {
public:
  NormalFence(const BSCurve &outer, const BSCurve &inner)
    : outer(outer), inner(inner) {
  }
  Point3D eval(double u) const override {
    VectorVector der;
    outer.eval(u, 1, der);
    auto q = inner.eval(u);
    Plane p;
    p.p = der[0];
    p.n = ((q - der[0]) ^ der[1]).normalize();
    return planeToSPoint(p);
  }
  Vector3D evalDerivative(double u) const override {
    // only used at u=0 and u=1
    if (u == 0.0)
      return (eval(epsilon) - eval(0.0)) / epsilon;
    return (eval(1.0) - eval(1.0 - epsilon)) / epsilon;
  }
private:
  BSCurve outer, inner;
};

Point3D readPoint(std::istream &is) {
  Point3D p;
  is >> p[0] >> p[1] >> p[2];
  return p;
}

BSCurve readCurve(std::istream &is) {
  size_t d, n_knots, n_points;
  is >> d >> n_points;
  n_knots = n_points + d + 1;

  DoubleVector knots(n_knots);
  for (size_t i = 0; i < n_knots; ++i)
    is >> knots[i];

  PointVector cpts(n_points);
  for (size_t i = 0; i < n_points; ++i)
    is >> cpts[i];

  return { d, knots, cpts };
}

auto readBoundaries(std::string filename) {
  std::ifstream f(filename.c_str());
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t n;
  f >> n;
  std::vector<std::shared_ptr<CurveType>> result;

  for (size_t i = 0; i < n; ++i) {
    auto outer = readCurve(f); outer.normalize();
    auto inner = readCurve(f); inner.normalize();
    result.push_back(std::make_shared<NormalFence>(outer, inner));
  }

  return result;
}


// Main code

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <model.dlp> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 15;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  auto boundaries = readBoundaries(argv[1]);
  C0Coons surface(boundaries);

  auto mesh = surface.eval(resolution);
  auto params = surface.parameters(resolution);
  auto n = params.size();
  std::vector<Vector3D> normals(n);

  for (size_t i = 0; i < n; ++i) {
    if (surface.onEdge(resolution, i)) {
      // On edge -> point and normal known a priori
      mesh[i] = Point3D(0,0,0);
      // TODO
    } else {
      // Not on edge -> normal from the mesh point, point from intersection
      auto p = spointToPlane(mesh[i]);
      auto p1 = spointToPlane(surface.eval(params[i] + Vector2D(epsilon, 0)));
      auto p2 = spointToPlane(surface.eval(params[i] + Vector2D(0, epsilon)));
      auto l1 = intersectPlanes(p, p1);
      auto l2 = intersectPlanes(p, p2);
      mesh[i] = intersectLines(l1, l2);
      normals[i] = p.n;
    }
  }

  auto sik = spointToPlane(mesh[0]);
  std::cout << "p: " << sik.p << std::endl;
  std::cout << "n: " << sik.n << std::endl;
  mesh.writeOBJ("/tmp/surface.obj");
}