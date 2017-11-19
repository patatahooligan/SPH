#pragma once

class Vec3f {
	public:
		float x, y, z;

		Vec3f(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}

		Vec3f(const Vec3f &other) = default;

		Vec3f& operator=(const Vec3f &other) = default;

		Vec3f& operator+=(const Vec3f &other) {
			x += other.x;
			y += other.y;
			z += other.z;
			return *this;
		}

		Vec3f& operator-=(const Vec3f &other) {
			x -= other.x;
			y -= other.y;
			z -= other.z;
			return *this;
		}

		Vec3f& operator*=(const float c) {
			x *= c;
			y *= c;
			z *= c;
			return *this;
		}

		Vec3f& operator/=(const float c) {
			x /= c;
			y /= c;
			z /= c;
			return *this;
		}

		float dot_product(const Vec3f &other) const {
			return x*other.x + y * other.y + z * other.z;
		}

		float length_squared() const {
			return x*x + y*y + z*z;
		}
		
		float length() const {
			return sqrt(length_squared());
		}
};

inline Vec3f operator+(Vec3f v1, const Vec3f &v2) {
	return v1 += v2;
}

inline Vec3f operator-(Vec3f v1, const Vec3f &v2) {
	return v1 -= v2;
}

inline Vec3f operator*(Vec3f v, const float c) {
	return v *= c;
}

inline Vec3f operator*(const float c, Vec3f v) {
	return v *= c;
}

inline Vec3f operator/(Vec3f v, const float c) {
	return v /= c;
}