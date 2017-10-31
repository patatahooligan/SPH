#include "stdafx.h"

#include <cmath>

#include "vec3f.h"


Vec3f::Vec3f(float setx, float sety, float setz) {
	x = setx;
	y = sety;
	z = setz;
}

Vec3f& Vec3f::operator+=(const Vec3f &other) {
	x += other.x;
	y += other.y;
	z += other.z;
	return *this;
}

Vec3f& Vec3f::operator-=(const Vec3f &other) {
	x -= other.x;
	y -= other.y;
	z -= other.z;
	return *this;
}

float Vec3f::dot_product(const Vec3f &other) const {
	return x*other.x + y * other.y + z * other.z;
}

float Vec3f::length_squared() const {
	return x*x + y*y + z*z;
}

float Vec3f::length() const {
	return sqrt(length_squared());
}


Vec3f operator+(Vec3f v1, const Vec3f &v2) {
	return v1 += v2;
}

Vec3f operator-(Vec3f v1, const Vec3f &v2) {
	return v1 -= v2;
}

Vec3f operator*(const Vec3f &v, const float c) {
	return Vec3f(
		v.x * c,
		v.y * c,
		v.z * c);
}

Vec3f operator*(const float c, const Vec3f &v) {
	return v * c;
}

Vec3f operator/(const Vec3f &v, const float c) {
	return Vec3f(
		v.x / c,
		v.y / c,
		v.z / c);
}