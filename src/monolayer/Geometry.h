#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
#include <cmath>
#include <iostream>

template <typename T> struct vec2d {
public:
	union {
		struct {
			T x, y;
		};
		T coord[2];
	};
	vec2d() {}
	vec2d(T x_, T y_) : x(x_), y(y_) {}
	vec2d(const vec2d<T>& other) : x(other.x), y(other.y) {}

	inline vec2d<T>& operator+=(const vec2d<T>& other) {
		x += other.x;
		y += other.y;
		return *this;
	}

	inline vec2d<T>& operator*=(const T scalar) {
		x *= scalar;
		y *= scalar;
		return *this;
	}

	inline vec2d<T>& operator/=(const T scalar) {
		x /= scalar;
		y /= scalar;
		return *this;
	}

	inline vec2d<T> operator+(const vec2d<T>& other) const {
		return vec2d<T>(this->x + other.x, this->y + other.y);
	}

	inline vec2d<T> operator-(const vec2d<T>& other) const {
		return vec2d<T>(this->x - other.x, this->y - other.y);
	}

	inline vec2d<T> operator*(const T scalar) const {
		return vec2d<T>(this->x * scalar, this->y * scalar);
	}

	inline vec2d<T> operator/(const T scalar) const {
		return vec2d<T>(this->x / scalar, this->y / scalar);
	}

	inline T operator*(const vec2d<T>& other) const {
		return this->x * other.x + this->y * other.y;
	}

	inline T norm2() const { return (*this) * (*this); }

	inline T norm() const { return std::sqrt(this->norm2()); }

	inline vec2d<T> normalized() const { return (*this) / this->norm(); }

	bool operator<(const vec2d<T>& rhs) const {
		if (x == rhs.x) {
			if (y == rhs.y) {
				return false;
			}
			else {
				return y < rhs.y;
			}
		}
		else {
			return x < rhs.x;
		}
	}

	bool operator==(const vec2d<T>& rhs) const {
		return y == rhs.y && x == rhs.x;
	}

	void print(std::ostream& stream) const {
		stream << "[" << x << "|" << y << "]";
	}
};

template <typename T>
std::ostream& operator<<(std::ostream& stream, const vec2d<T>& entry) {
	entry.print(stream);
	return stream;
}



template <typename T> struct vec3d {
public:
	//T x, y, z;
	union {
		struct {
			T x, y, z;
		};
		T coord[3];
	};
	vec3d() {}
	vec3d(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}
	vec3d(const vec3d<T>& other) : x(other.x), y(other.y), z(other.z) {}

	inline vec3d<T>& operator+=(const vec3d<T>& other) {
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	inline vec3d<T>& operator*=(const T scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	inline vec3d<T>& operator/=(const T scalar) {
		x /= scalar;
		y /= scalar;
		z /= scalar;
		return *this;
	}

	inline vec3d<T> operator+(const vec3d<T>& other) const {
		return vec3d<T>(this->x + other.x, this->y + other.y,
			this->z + other.z);
	}

	inline vec3d<T> operator-(const vec3d<T>& other) const {
		return vec3d<T>(this->x - other.x, this->y - other.y,
			this->z - other.z);
	}

	inline vec3d<T> operator*(const T scalar) const {
		return vec3d<T>(this->x * scalar, this->y * scalar, this->z * scalar);
	}

	inline vec3d<T> operator/(const T scalar) const {
		return vec3d<T>(this->x / scalar, this->y / scalar, this->z / scalar);
	}

	inline T operator*(const vec3d<T>& other) const {
		return this->x * other.x + this->y * other.y + this->z * other.z;
	}

	inline vec3d<T> cross(const vec3d<T>& other) const {
		return vec3d<T>(this->y * other.z - this->z * other.y,
			this->z * other.x - this->x * other.z,
			this->x * other.y - this->y * other.x);
	}

	inline T norm2() const { return (*this) * (*this); }

	inline T norm() const { return std::sqrt(this->norm2()); }

	inline vec3d<T> normalized() const { return (*this) / this->norm(); }

	bool operator<(const vec3d<T>& rhs) const {
		if (x == rhs.x) {
			if (y == rhs.y) {
				if (z == rhs.z) {
					return false;
				}
				else {
					return z < rhs.z;
				}
			}
			else {
				return y < rhs.y;
			}
		}
		else {
			return x < rhs.x;
		}
	}

	bool operator==(const vec3d<T>& rhs) const {
		return z == rhs.z && y == rhs.y && x == rhs.x;
	}

	void print(std::ostream& stream) const {
		stream << "[" << x << "|" << y << "|" << z << "]";
	}
};

template <typename T>
std::ostream& operator<<(std::ostream& stream, const vec3d<T>& entry) {
	entry.print(stream);
	return stream;
}
#endif
