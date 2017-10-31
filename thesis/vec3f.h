#pragma once

class Vec3f {
	public:
		float x, y, z;

		Vec3f(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}

		Vec3f(const Vec3f &other) = default;

		Vec3f& operator=(const Vec3f &other) = default;

		Vec3f& operator+=(const Vec3f &other);

		Vec3f& operator-=(const Vec3f &other);

		float dot_product(const Vec3f &other) const;

		float length_squared() const;
		
		float length() const;
};

Vec3f operator+(Vec3f v1, const Vec3f &v2);

Vec3f operator-(Vec3f v1, const Vec3f &v2);

Vec3f operator*(const Vec3f &v, const float c);

Vec3f operator*(const float c, const Vec3f &v);

Vec3f operator/(const Vec3f &v, const float c);