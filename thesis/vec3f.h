#pragma once

// Two declarations to resolve the circular dependency issue between class methods and free functions
class Vec3f;
inline Vec3f operator/(Vec3f v, const float c);

class Vec3f {
	public:
		float x, y, z;

		constexpr Vec3f(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}

		constexpr Vec3f(const Vec3f &other) = default;

		Vec3f& operator=(const Vec3f &other) = default;

		bool operator==(const Vec3f &other) const {
			return x == other.x && y == other.y && z == other.z;
		}

		bool operator!=(const Vec3f &other) const {
			return !(*this == other);
		}

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

		// Calling this with out of range i is undefined behavior
		// Ignore the warning that not all control paths return a value
		#pragma warning(suppress: 4715)
		const float& operator[](size_t i) const {
			switch (i) {
			case 0:
				return x;
			case 1:
				return y;
			case 2:
				return z;
			}
		}

		float& operator[](size_t i) {
			return const_cast<float&>(const_cast<const Vec3f&>(*this)[i]);
		}

		float dot_product(const Vec3f &other) const {
			return x*other.x + y * other.y + z * other.z;
		}

		float length_squared() const {
			return x*x + y*y + z*z;
		}
		
		float length() const {
			return std::sqrt(length_squared());
		}

		Vec3f unit_vector() const {
			return *this / this->length();
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

inline float dot_product(const Vec3f& v1, const Vec3f& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline bool isfinite(Vec3f v) {
	return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}