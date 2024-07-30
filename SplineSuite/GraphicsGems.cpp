/* GGVecLib.c */
/*
2d and 3d Vector C Library
by Andrew Glassner
from "Graphics Gems", Academic Press, 1990
*/

#include "GraphicsGems.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/******************/
/*   2d Library   */
/******************/

/* returns squared length of input vector */
double V2SquaredLength(Vector2 *a) { return ((a->x * a->x) + (a->y * a->y)); }

/* returns length of input vector */
double V2Length(Vector2 *a) { return sqrt(V2SquaredLength(a)); }

/* negates the input vector and returns it */
Vector2 *V2Negate(Vector2 *v) {
  v->x = -v->x;
  v->y = -v->y;
  return v;
}

/* normalizes the input vector and returns it */
Vector2 *V2Normalize(Vector2 *v) {
  double len = V2Length(v);
  if (len != 0.0) {
    v->x /= len;
    v->y /= len;
  }
  return v;
}

/* scales the input vector to the new length and returns it */
Vector2 *V2Scale(Vector2 *v, double newlen) {
  double len = V2Length(v);
  if (len != 0.0) {
    v->x *= newlen / len;
    v->y *= newlen / len;
  }
  return v;
}

/* return vector sum c = a+b */
Vector2 *V2Add(Vector2 *a, Vector2 *b, Vector2 *c) {
  c->x = a->x + b->x;
  c->y = a->y + b->y;
  return c;
}

/* return vector difference c = a-b */
Vector2 *V2Sub(Vector2 *a, Vector2 *b, Vector2 *c) {
  c->x = a->x - b->x;
  c->y = a->y - b->y;
  return c;
}

/* return the dot product of vectors a and b */
double V2Dot(Vector2 *a, Vector2 *b) { return (a->x * b->x) + (a->y * b->y); }

/* linearly interpolate between vectors by an amount alpha */
/* and return the resulting vector. */
/* When alpha=0, result=lo.  When alpha=1, result=hi. */
Vector2 *V2Lerp(Vector2 *lo, Vector2 *hi, double alpha, Vector2 *result) {
  result->x = LERP(alpha, lo->x, hi->x);
  result->y = LERP(alpha, lo->y, hi->y);
  return result;
}

/* make a linear combination of two vectors and return the result. */
/* result = (a * ascl) + (b * bscl) */
Vector2 *V2Combine(Vector2 *a, Vector2 *b, Vector2 *result, double ascl,
                   double bscl) {
  result->x = (ascl * a->x) + (bscl * b->x);
  result->y = (ascl * a->y) + (bscl * b->y);
  return result;
}

/* multiply two vectors together component-wise */
Vector2 *V2Mul(Vector2 *a, Vector2 *b, Vector2 *result) {
  result->x = a->x * b->x;
  result->y = a->y * b->y;
  return result;
}

/* return the distance between two points */
double V2DistanceBetween2Points(Point2 *a, Point2 *b) {
  double dx = a->x - b->x;
  double dy = a->y - b->y;
  return sqrt((dx * dx) + (dy * dy));
}

/* return the vector perpendicular to the input vector a */
Vector2 *V2MakePerpendicular(Vector2 *a, Vector2 *ap) {
  ap->x = -a->y;
  ap->y = a->x;
  return ap;
}

/* create, initialize, and return a new vector */
Vector2 *V2New(double x, double y) {
  Vector2 *v = NEWTYPE(Vector2);
  v->x = x;
  v->y = y;
  return v;
}

/* create, initialize, and return a duplicate vector */
Vector2 *V2Duplicate(Vector2 *a) {
  Vector2 *v = NEWTYPE(Vector2);
  v->x = a->x;
  v->y = a->y;
  return v;
}

/* multiply a point by a projective matrix and return the transformed point */
Point2 *V2MulPointByProjMatrix(Point2 *pin, Matrix3 *m, Point2 *pout) {
  double w;
  pout->x = (pin->x * m->element[0][0]) + (pin->y * m->element[1][0]) +
            m->element[2][0];
  pout->y = (pin->x * m->element[0][1]) + (pin->y * m->element[1][1]) +
            m->element[2][1];
  w = (pin->x * m->element[0][2]) + (pin->y * m->element[1][2]) +
      m->element[2][2];
  if (w != 0.0) {
    pout->x /= w;
    pout->y /= w;
  }
  return pout;
}

/* multiply together matrices c = ab */
/* note that c must not point to either of the input matrices */
Matrix3 *V2MatMul(Matrix3 *a, Matrix3 *b, Matrix3 *c) {
  int i, j, k;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      c->element[i][j] = 0;
      for (k = 0; k < 3; k++) {
        c->element[i][j] += a->element[i][k] * b->element[k][j];
      }
    }
  }
  return c;
}

/* transpose matrix a, return b */
Matrix3 *TransposeMatrix3(Matrix3 *a, Matrix3 *b) {
  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      b->element[i][j] = a->element[j][i];
    }
  }
  return b;
}

/* returns squared length of input vector */
double V3SquaredLength(Vector3 *a) {
  return (a->x * a->x) + (a->y * a->y) + (a->z * a->z);
}

/* returns length of input vector */
double V3Length(Vector3 *a) { return sqrt(V3SquaredLength(a)); }

/* negates the input vector and returns it */
Vector3 *V3Negate(Vector3 *v) {
  v->x = -v->x;
  v->y = -v->y;
  v->z = -v->z;
  return v;
}

/* normalizes the input vector and returns it */
Vector3 *V3Normalize(Vector3 *v) {
  double len = V3Length(v);
  if (len != 0.0) {
    v->x /= len;
    v->y /= len;
    v->z /= len;
  }
  return v;
}

/* scales the input vector to the new length and returns it */
Vector3 *V3Scale(Vector3 *v, double newlen) {
  double len = V3Length(v);
  if (len != 0.0) {
    v->x *= newlen / len;
    v->y *= newlen / len;
    v->z *= newlen / len;
  }
  return v;
}

/* return vector sum c = a+b */
Vector3 *V3Add(Vector3 *a, Vector3 *b, Vector3 *c) {
  c->x = a->x + b->x;
  c->y = a->y + b->y;
  c->z = a->z + b->z;
  return c;
}

/* return vector difference c = a-b */
Vector3 *V3Sub(Vector3 *a, Vector3 *b, Vector3 *c) {
  c->x = a->x - b->x;
  c->y = a->y - b->y;
  c->z = a->z - b->z;
  return c;
}

/* return the dot product of vectors a and b */
double V3Dot(Vector3 *a, Vector3 *b) {
  return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

/* linearly interpolate between vectors by an amount alpha */
/* and return the resulting vector. */
/* When alpha=0, result=lo.  When alpha=1, result=hi. */
Vector3 *V3Lerp(Vector3 *lo, Vector3 *hi, double alpha, Vector3 *result) {
  result->x = LERP(alpha, lo->x, hi->x);
  result->y = LERP(alpha, lo->y, hi->y);
  result->z = LERP(alpha, lo->z, hi->z);
  return result;
}

/* make a linear combination of two vectors and return the result. */
/* result = (a * ascl) + (b * bscl) */
Vector3 *V3Combine(Vector3 *a, Vector3 *b, Vector3 *result, double ascl,
                   double bscl) {
  result->x = (ascl * a->x) + (bscl * b->x);
  result->y = (ascl * a->y) + (bscl * b->y);
  result->z = (ascl * a->z) + (bscl * b->z);
  return result;
}

/* multiply two vectors together component-wise and return the result */
Vector3 *V3Mul(Vector3 *a, Vector3 *b, Vector3 *result) {
  result->x = a->x * b->x;
  result->y = a->y * b->y;
  result->z = a->z * b->z;
  return result;
}

/* return the distance between two points */
double V3DistanceBetween2Points(Point3 *a, Point3 *b) {
  double dx = a->x - b->x;
  double dy = a->y - b->y;
  double dz = a->z - b->z;
  return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

/* return the cross product c = a cross b */
Vector3 *V3Cross(Vector3 *a, Vector3 *b, Vector3 *c) {
  c->x = (a->y * b->z) - (a->z * b->y);
  c->y = (a->z * b->x) - (a->x * b->z);
  c->z = (a->x * b->y) - (a->y * b->x);
  return c;
}

/* create, initialize, and return a new vector */
Vector3 *V3New(double x, double y, double z) {
  Vector3 *v = NEWTYPE(Vector3);
  v->x = x;
  v->y = y;
  v->z = z;
  return v;
}

/* create, initialize, and return a duplicate vector */
Vector3 *V3Duplicate(Vector3 *a) {
  Vector3 *v = NEWTYPE(Vector3);
  v->x = a->x;
  v->y = a->y;
  v->z = a->z;
  return v;
}

/* multiply a point by a matrix and return the transformed point */
Point3 *V3MulPointByMatrix(Point3 *pin, Matrix3 *m, Point3 *pout) {
  pout->x = (pin->x * m->element[0][0]) + (pin->y * m->element[1][0]) +
            (pin->z * m->element[2][0]);
  pout->y = (pin->x * m->element[0][1]) + (pin->y * m->element[1][1]) +
            (pin->z * m->element[2][1]);
  pout->z = (pin->x * m->element[0][2]) + (pin->y * m->element[1][2]) +
            (pin->z * m->element[2][2]);
  return pout;
}

/* multiply a point by a projective matrix and return the transformed point */
Point3 *V3MulPointByProjMatrix(Point3 *pin, Matrix4 *m, Point3 *pout) {
  double w;
  pout->x = (pin->x * m->element[0][0]) + (pin->y * m->element[1][0]) +
            (pin->z * m->element[2][0]) + m->element[3][0];
  pout->y = (pin->x * m->element[0][1]) + (pin->y * m->element[1][1]) +
            (pin->z * m->element[2][1]) + m->element[3][1];
  pout->z = (pin->x * m->element[0][2]) + (pin->y * m->element[1][2]) +
            (pin->z * m->element[2][2]) + m->element[3][2];
  w = (pin->x * m->element[0][3]) + (pin->y * m->element[1][3]) +
      (pin->z * m->element[2][3]) + m->element[3][3];
  if (w != 0.0) {
    pout->x /= w;
    pout->y /= w;
    pout->z /= w;
  }
  return pout;
}

/* multiply together matrices c = ab */
/* note that c must not point to either of the input matrices */
Matrix4 *V3MatMul(Matrix4 *a, Matrix4 *b, Matrix4 *c) {
  int i, j, k;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      c->element[i][j] = 0;
      for (k = 0; k < 4; k++) {
        c->element[i][j] += a->element[i][k] * b->element[k][j];
      }
    }
  }
  return c;
}

/* binary greatest common divisor by Silver and Terzian.  See Knuth */
/* both inputs must be >= 0 */
int gcd(int u, int v) {
  int t, f;
  if (u < 0 || v < 0)
    return 1; /* error if u<0 or v<0 */
  f = 1;
  while ((u % 2 == 0) && (v % 2 == 0)) {
    u >>= 1;
    v >>= 1;
    f *= 2;
  }
  if (u & 1) {
    t = -v;
    goto B4;
  } else {
    t = u;
  }
B3:
  if (t > 0) {
    t >>= 1;
  } else {
    t = -((-t) >> 1);
  }
B4:
  if (t % 2 == 0)
    goto B3;

  if (t > 0) {
    u = t;
  } else {
    v = -t;
  }
  if ((t = u - v) != 0)
    goto B3;
  return u * f;
}

/* return roots of ax^2+bx+c */
/* stable algebra derived from Numerical Recipes by Press et al.*/
int quadraticRoots(double a, double b, double c, double *roots) {
  double d, q;
  int count = 0;
  d = (b * b) - (4 * a * c);
  if (d < 0.0) {
    *roots = *(roots + 1) = 0.0;
    return 0;
  }
  q = -0.5 * (b + (SGN(b) * sqrt(d)));
  if (a != 0.0) {
    *roots++ = q / a;
    count++;
  }
  if (q != 0.0) {
    *roots++ = c / q;
    count++;
  }
  return count;
}

/* generic 1d regula-falsi step.  f is function to evaluate */
/* interval known to contain root is given in left, right */
/* returns new estimate */
double RegulaFalsi(double (*f)(double), double left, double right) {
  double d = (*f)(right) - (*f)(left);
  if (d != 0.0) {
    return (right - (*f)(right) * (right - left) / d);
  }
  return (left + right) / 2.0;
}

/* generic 1d Newton-Raphson step. f is function, df is derivative */
/* x is current best guess for root location. Returns new estimate */
double NewtonRaphson(double (*f)(double), double (*df)(double), double x) {
  double d = (*df)(x);
  if (d != 0.0) {
    return (x - ((*f)(x) / d));
  }
  return x - 1.0;
}

/* hybrid 1d Newton-Raphson/Regula Falsi root finder. */
/* input function f and its derivative df, an interval */
/* left, right known to contain the root, and an error tolerance */
/* Based on Blinn */
double findroot(double left, double right, double tolerance,
                double (*f)(double), double (*df)(double)) {
  double newx = left;
  while (ABS((*f)(newx)) > tolerance) {
    newx = NewtonRaphson(f, df, newx);
    if (newx < left || newx > right) {
      newx = RegulaFalsi(f, left, right);
    }
    if ((*f)(newx) * (*f)(left) <= 0.0) {
      right = newx;
    } else {
      left = newx;
    }
  }
  return newx;
}

/*
Solving the Nearest Point-on-Curve Problem
and
A Bezier Curve-Based Root-Finder
by Philip J. Schneider
from "Graphics Gems", Academic Press, 1990
*/

/*	point_on_curve.c	*/

#define TESTMODE

/*
 *  Forward declarations
 */
// Point2  NearestPointOnCurve();
// static	int	FindRoots();
// static	Point2	*ConvertToBezierForm();
// static	double	ComputeXIntercept();
// static	int	ControlPolygonFlatEnough();
// static	int	CrossingCount();
// static	Point2	Bezier();
// static	Vector2	V2ScaleII();

// #ifdef TESTMODE
// /*
//  *  main :
//  *	Given a cubic Bezier curve (i.e., its control points), and some
//  *	arbitrary point in the plane, find the point on the curve
//  *	closest to that arbitrary point.
//  */
// int main()
// {

//  static Point2 bezCurve[4] = {	/*  A cubic Bezier curve	*/
// 	{ 0.0, 0.0 },
// 	{ 1.0, 2.0 },
// 	{ 3.0, 3.0 },
// 	{ 4.0, 2.0 },
//     };
//     static Point2 arbPoint = { 3.5, 2.0 }; /*Some arbitrary point*/
//     Point2	pointOnCurve;		 /*  Nearest point on the curve */

//     /*  Find the closest point */
//     pointOnCurve = NearestPointOnCurve(arbPoint, bezCurve);
//     printf("pointOnCurve : (%4.4f, %4.4f)\n", pointOnCurve.x,
// 		pointOnCurve.y);
// }
// #endif /* TESTMODE */

/*
 *  Bezier :
 *	Evaluate a Bezier curve at a particular parameter value
 *      Fill in control points for resulting sub-curves if "Left" and
 *	"Right" are non-null.
 *
 */
Point2 Bezier(Point2 *V, int degree, double t, Point2 *Left, Point2 *Right)

{
  int i, j; /* Index variables	*/
  Point2 Vtemp[W_DEGREE + 1][W_DEGREE + 1];

  /* Copy control points	*/
  for (j = 0; j <= degree; j++) {
    Vtemp[0][j] = V[j];
  }

  /* Triangle computation	*/
  for (i = 1; i <= degree; i++) {
    for (j = 0; j <= degree - i; j++) {
      Vtemp[i][j].x = (1.0 - t) * Vtemp[i - 1][j].x + t * Vtemp[i - 1][j + 1].x;
      Vtemp[i][j].y = (1.0 - t) * Vtemp[i - 1][j].y + t * Vtemp[i - 1][j + 1].y;
    }
  }

  if (Left != NULL) {
    for (j = 0; j <= degree; j++) {
      Left[j] = Vtemp[j][0];
    }
  }
  if (Right != NULL) {
    for (j = 0; j <= degree; j++) {
      Right[j] = Vtemp[degree - j][j];
    }
  }

  return (Vtemp[degree][0]);
}

/* scales the input vector by the scalar s and returns the result */
Vector2 V2ScaleII(Vector2 *v, double s) {
  Vector2 result;

  result.x = v->x * s;
  result.y = v->y * s;
  return result;
}

/*
 *  ComputeXIntercept :
 *	Compute intersection of chord from first control point to last
 *  	with 0-axis.
 *
 */
/* NOTE: "T" and "Y" do not have to be computed, and there are many useless
 * operations in the following (e.g. "0.0 - 0.0").
 */
double ComputeXIntercept(Point2 *V, int degree)

{
  double XLK, YLK, XNM, YNM, XMK, YMK;
  double det, detInv;
  double S;
  double X;

  XLK = 1.0 - 0.0;
  YLK = 0.0 - 0.0;
  XNM = V[degree].x - V[0].x;
  YNM = V[degree].y - V[0].y;
  XMK = V[0].x - 0.0;
  YMK = V[0].y - 0.0;

  det = XNM * YLK - YNM * XLK;
  detInv = 1.0 / det;

  S = (XNM * YMK - YNM * XMK) * detInv;
  /*  T = (XLK*YMK - YLK*XMK) * detInv; */

  X = 0.0 + XLK * S;
  /*  Y = 0.0 + YLK * S; */

  return X;
}

/*
 *  ControlPolygonFlatEnough :
 *	Check if the control polygon of a Bezier curve is flat enough
 *	for recursive subdivision to bottom out.
 *
 *  Corrections by James Walker, jw@jwwalker.com, as follows:

There seem to be errors in the ControlPolygonFlatEnough function in the
Graphics Gems book and the repository (NearestPoint.c). This function
is briefly described on p. 413 of the text, and appears on pages 793-794.
I see two main problems with it.

The idea is to find an upper bound for the error of approximating the x
intercept of the Bezier curve by the x intercept of the line through the
first and last control points. It is claimed on p. 413 that this error is
bounded by half of the difference between the intercepts of the bounding
box. I don't see why that should be true. The line joining the first and
last control points can be on one side of the bounding box, and the actual
curve can be near the opposite side, so the bound should be the difference
of the bounding box intercepts, not half of it.

Second, we come to the implementation. The values distance[i] computed in
the first loop are not actual distances, but squares of distances. I
realize that minimizing or maximizing the squares is equivalent to
minimizing or maximizing the distances.  But when the code claims that
one of the sides of the bounding box has equation
a * x + b * y + c + max_distance_above, where max_distance_above is one of
those squared distances, that makes no sense to me.

I have appended my version of the function. If you apply my code to the
cubic Bezier curve used to test NearestPoint.c,

 static Point2 bezCurve[4] = {    /  A cubic Bezier curve    /
    { 0.0, 0.0 },
    { 1.0, 2.0 },
    { 3.0, 3.0 },
    { 4.0, 2.0 },
    };

my code computes left_intercept = -3.0 and right_intercept = 0.0, which you
can verify by sketching a graph. The original code computes
left_intercept = 0.0 and right_intercept = 0.9.

 */

/* static int ControlPolygonFlatEnough( const Point2* V, int degree ) */
int ControlPolygonFlatEnough(Point2 *V, int degree)

{
  int i; /* Index variable        */
  double value;
  double max_distance_above;
  double max_distance_below;
  double error; /* Precision of root        */
  double intercept_1, intercept_2, left_intercept, right_intercept;
  double a, b, c; /* Coefficients of implicit    */
                  /* eqn for line from V[0]-V[deg]*/
  double det, dInv;
  double a1, b1, c1, a2, b2, c2;

  /* Derive the implicit equation for line connecting first and last control
   * points */
  a = V[0].y - V[degree].y;
  b = V[degree].x - V[0].x;
  c = V[0].x * V[degree].y - V[degree].x * V[0].y;

  max_distance_above = max_distance_below = 0.0;

  for (i = 1; i < degree; i++) {
    value = a * V[i].x + b * V[i].y + c;

    if (value > max_distance_above) {
      max_distance_above = value;
    } else if (value < max_distance_below) {
      max_distance_below = value;
    }
  }

  /*  Implicit equation for zero line */
  a1 = 0.0;
  b1 = 1.0;
  c1 = 0.0;

  /*  Implicit equation for "above" line */
  a2 = a;
  b2 = b;
  c2 = c - max_distance_above;

  det = a1 * b2 - a2 * b1;
  dInv = 1.0 / det;

  intercept_1 = (b1 * c2 - b2 * c1) * dInv;

  /*  Implicit equation for "below" line */
  a2 = a;
  b2 = b;
  c2 = c - max_distance_below;

  det = a1 * b2 - a2 * b1;
  dInv = 1.0 / det;

  intercept_2 = (b1 * c2 - b2 * c1) * dInv;

  /* Compute intercepts of bounding box    */
  left_intercept = MIN(intercept_1, intercept_2);
  right_intercept = MAX(intercept_1, intercept_2);

  error = right_intercept - left_intercept;

  return (error < EPSILON) ? 1 : 0;
}

/*
 * CrossingCount :
 *	Count the number of times a Bezier control polygon
 *	crosses the 0-axis. This number is >= the number of roots.
 *
 */
int CrossingCount(Point2 *V, int degree)

{
  int i;
  int n_crossings = 0; /*  Number of zero-crossings	*/
  int sign, old_sign;  /*  Sign of coefficients	*/

  sign = old_sign = SGN(V[0].y);
  for (i = 1; i <= degree; i++) {
    sign = SGN(V[i].y);
    if (sign != old_sign)
      n_crossings++;
    old_sign = sign;
  }
  return n_crossings;
}

/*
 *  FindRoots :
 *	Given a 5th-degree equation in Bernstein-Bezier form, find
 *	all of the roots in the interval [0, 1].  Return the number
 *	of roots found.
 */
int FindRoots(Point2 *w, int degree, double *t, int depth) {
  int i;
  Point2 Left[W_DEGREE + 1],   /* New left and right 		*/
      Right[W_DEGREE + 1];     /* control polygons		*/
  int left_count,              /* Solution count from		*/
      right_count;             /* children			*/
  double left_t[W_DEGREE + 1], /* Solutions from kids		*/
      right_t[W_DEGREE + 1];

  switch (CrossingCount(w, degree)) {
  case 0: { /* No solutions here	*/
    return 0;
  }
  case 1: { /* Unique solution	*/
    /* Stop recursion when the tree is deep enough	*/
    /* if deep enough, return 1 solution at midpoint 	*/
    if (depth >= MAXDEPTH) {
      t[0] = (w[0].x + w[W_DEGREE].x) / 2.0;
      return 1;
    }
    if (ControlPolygonFlatEnough(w, degree)) {
      t[0] = ComputeXIntercept(w, degree);
      return 1;
    }
    break;
  }
  }

  /* Otherwise, solve recursively after	*/
  /* subdividing control polygon		*/
  Bezier(w, degree, 0.5, Left, Right);
  left_count = FindRoots(Left, degree, left_t, depth + 1);
  right_count = FindRoots(Right, degree, right_t, depth + 1);

  /* Gather solutions together	*/
  for (i = 0; i < left_count; i++) {
    t[i] = left_t[i];
  }
  for (i = 0; i < right_count; i++) {
    t[i + left_count] = right_t[i];
  }

  /* Send back total number of solutions	*/
  return (left_count + right_count);
}

/*
 *  ConvertToBezierForm :
 *		Given a point and a Bezier curve, generate a 5th-degree
 *		Bezier-format equation whose solution finds the point on the
 *      curve nearest the user-defined point.
 */

Point2 *ConvertToBezierForm(Point2 P, Point2 *V) {
  int i, j, k, m, n, ub, lb;
  int row, column;       /* Table indices		*/
  Vector2 c[DEGREE + 1]; /* V(i)'s - P			*/
  Vector2 d[DEGREE];     /* V(i+1) - V(i)		*/
  Point2 *w;             /* Ctl pts of 5th-degree curve  */
  double cdTable[3][4];  /* Dot product of c, d		*/
  static double z[3][4] = {
      /* Precomputed "z" for cubics	*/
      {1.0, 0.6, 0.3, 0.1},
      {0.4, 0.6, 0.6, 0.4},
      {0.1, 0.3, 0.6, 1.0},
  };

  /*Determine the c's -- these are vectors created by subtracting*/
  /* point P from each of the control points				*/
  for (i = 0; i <= DEGREE; i++) {
    V2Sub(&V[i], &P, &c[i]);
  }
  /* Determine the d's -- these are vectors created by subtracting*/
  /* each control point from the next					*/
  for (i = 0; i <= DEGREE - 1; i++) {
    d[i] = V2ScaleII(V2Sub(&V[i + 1], &V[i], &d[i]), 3.0);
  }

  /* Create the c,d table -- this is a table of dot products of the */
  /* c's and d's							*/
  for (row = 0; row <= DEGREE - 1; row++) {
    for (column = 0; column <= DEGREE; column++) {
      cdTable[row][column] = V2Dot(&d[row], &c[column]);
    }
  }

  /* Now, apply the z's to the dot products, on the skew diagonal*/
  /* Also, set up the x-values, making these "points"		*/
  w = (Point2 *)malloc((unsigned)(W_DEGREE + 1) * sizeof(Point2));
  for (i = 0; i <= W_DEGREE; i++) {
    w[i].y = 0.0;
    w[i].x = (double)(i) / W_DEGREE;
  }

  n = DEGREE;
  m = DEGREE - 1;
  for (k = 0; k <= n + m; k++) {
    lb = MAX(0, k - m);
    ub = MIN(k, n);
    for (i = lb; i <= ub; i++) {
      j = k - i;
      w[i + j].y += cdTable[j][i] * z[j][i];
    }
  }

  return (w);
}

/*
 *  NearestPointOnCurve :
 *  	Compute the parameter value of the point on a Bezier
 *		curve segment closest to some arbtitrary, user-input point.
 *		Return the point on the curve at that parameter value.
 *
 */
Point2 NearestPointOnCurve(Point2 P, Point2 *V) {
  Point2 *w;                    /* Ctl pts for 5th-degree eqn	*/
  double t_candidate[W_DEGREE]; /* Possible roots		*/
  int n_solutions;              /* Number of roots found	*/
  double t;                     /* Parameter value of closest pt*/

  /*  Convert problem to 5th-degree Bezier form	*/
  w = ConvertToBezierForm(P, V);

  /* Find all possible roots of 5th-degree equation */
  n_solutions = FindRoots(w, W_DEGREE, t_candidate, 0);
  free((char *)w);

  /* Compare distances of P to all candidates, and to t=0, and t=1 */
  {
    double dist, new_dist;
    Point2 p;
    Vector2 v;
    int i;

    /* Check distance to beginning of curve, where t = 0	*/
    dist = V2SquaredLength(V2Sub(&P, &V[0], &v));
    t = 0.0;

    /* Find distances for candidate points	*/
    for (i = 0; i < n_solutions; i++) {
      p = Bezier(V, DEGREE, t_candidate[i], (Point2 *)NULL, (Point2 *)NULL);
      new_dist = V2SquaredLength(V2Sub(&P, &p, &v));
      if (new_dist < dist) {
        dist = new_dist;
        t = t_candidate[i];
      }
    }

    /* Finally, look at distance to end point, where t = 1.0 */
    new_dist = V2SquaredLength(V2Sub(&P, &V[DEGREE], &v));
    if (new_dist < dist) {
      dist = new_dist;
      t = 1.0;
    }
  }

  /*  Return the point on the curve at parameter value t */
  //     printf("t : %4.12f\n", t);
  return (Bezier(V, DEGREE, t, (Point2 *)NULL, (Point2 *)NULL));
}
