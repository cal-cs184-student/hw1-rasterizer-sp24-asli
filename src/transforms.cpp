#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

#include <cmath>
# define M_PI           3.14159265358979323846  /* pi */

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float dx, float dy) {
	// Part 3: Fill this in.

	return Matrix3x3(1, 0, dx, 0, 1, dy, 0, 0, 1);
}

Matrix3x3 scale(float sx, float sy) {
	// Elements of the matrix are filled in row-major order.
	// 1st row: sx 0 0
	// 2nd row: 0 sy 0
	// 3rd row: 0 0 1
	return Matrix3x3(sx, 0, 0,
		0, sy, 0,
		0, 0, 1);
}

// The input argument is in degrees counterclockwise
	Matrix3x3 rotate(float deg) {
		// Convert degrees to radians
		float rad = deg * (M_PI / 180.0);
		return Matrix3x3(cos(rad), -sin(rad), 0,
			sin(rad), cos(rad), 0,
			0, 0, 1);
	}
}
