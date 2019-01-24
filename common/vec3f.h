#pragma once

// Two declarations to resolve the circular dependency issue between class methods and free functions
class Vec3f;
constexpr inline Vec3f operator/(Vec3f v, const float c);

class Vec3f {
	public:
		float x, y, z;

		constexpr Vec3f(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}

		constexpr Vec3f(const Vec3f &other) = default;

		static constexpr Vec3f x_unit() { return { 1.0f, 0.0f, 0.0f }; }
		static constexpr Vec3f y_unit() { return { 0.0f, 1.0f, 0.0f }; }
		static constexpr Vec3f z_unit() { return { 0.0f, 0.0f, 1.0f }; }

		constexpr Vec3f& operator=(const Vec3f &other) = default;

		constexpr bool operator==(const Vec3f &other) const {
			return x == other.x && y == other.y && z == other.z;
		}

		constexpr bool operator!=(const Vec3f &other) const {
			return !(*this == other);
		}

		constexpr Vec3f& operator+=(const Vec3f &other) {
			x += other.x;
			y += other.y;
			z += other.z;
			return *this;
		}

		constexpr Vec3f& operator-=(const Vec3f &other) {
			x -= other.x;
			y -= other.y;
			z -= other.z;
			return *this;
		}

		constexpr Vec3f& operator*=(const float c) {
			x *= c;
			y *= c;
			z *= c;
			return *this;
		}

		constexpr Vec3f& operator/=(const float c) {
			x /= c;
			y /= c;
			z /= c;
			return *this;
		}

		// Calling this with out of range i is undefined behavior
		// Ignore the warning that not all control paths return a value
		const float& operator[](size_t i) const {
			switch (i) {
			case 0:
				return x;
			case 1:
				return y;
			case 2:
				return z;
			}
			__builtin_unreachable();
		}

		float& operator[](size_t i) {
			return const_cast<float&>(const_cast<const Vec3f&>(*this)[i]);
		}

		constexpr float dot_product(const Vec3f &other) const {
			return x*other.x + y * other.y + z * other.z;
		}

		constexpr float length_squared() const {
			return x*x + y*y + z*z;
		}
		
		float length() const {
			return std::sqrt(length_squared());
		}

		Vec3f unit_vector() const {
			return *this / this->length();
		}

		Vec3f rotate_x(const float angle_in_rad) const {
			return {
				x,
				y * std::cos(angle_in_rad) + z * std::sin(angle_in_rad),
				- y * std::sin(angle_in_rad) + z * std::cos(angle_in_rad) };
		}

		Vec3f rotate_y(const float angle_in_rad) const {
			return {
				x * std::cos(angle_in_rad) - z * std::sin(angle_in_rad),
				y,
				x * std::sin(angle_in_rad) + z * std::cos(angle_in_rad) };
		}

		Vec3f rotate_z(const float angle_in_rad) const {
			return {
				x * std::cos(angle_in_rad) + y * std::sin(angle_in_rad),
				-x * std::sin(angle_in_rad) + y * std::cos(angle_in_rad),
				z };
		}
};

constexpr inline Vec3f operator+(Vec3f v1, const Vec3f &v2) {
	return v1 += v2;
}

constexpr inline Vec3f operator-(Vec3f v1, const Vec3f &v2) {
	return v1 -= v2;
}

constexpr inline Vec3f operator-(const Vec3f &v) {
	return Vec3f{ 0.0f, 0.0f, 0.0f } - v;
}

constexpr inline Vec3f operator*(Vec3f v, const float c) {
	return v *= c;
}

constexpr inline Vec3f operator*(const float c, Vec3f v) {
	return v *= c;
}

constexpr inline Vec3f operator/(Vec3f v, const float c) {
	return v /= c;
}

constexpr inline float dot_product(const Vec3f& v1, const Vec3f& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline bool isfinite(Vec3f v) {
	return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}