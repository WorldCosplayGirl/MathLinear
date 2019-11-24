//===========================================================================
// Copyright (C)2019  Berk Atabek ( atabek dot berk at hotmail dot com )
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  File	    : MathLinear.h      
//  Desc	    : Small size linear math library for 3D graphics.
//
//                [Classes]
//                vec2d    : 2D vector data structure
//				  vec3d    : 3D vector data structure
//				  vec4d    : 4D vector data structure 
//                mat3x3   : 3x3 Matrix
//                mat4x4   : 4x4 Matrix
//                ray      : Parametric line segment
//                plane    : 3D plane structure
//				  rectangle: 2D bounding rectangle 
//				  AABB     : Axis-Aligned-Bounding-Box
//                frustum  : view frustum for visibility testing 
//				  colorRGB : RGB triplet color structure
//
//  Note        : All matrix transformations are right handed and row major.
//                In order to use them with OpenGL pipeline users need to use
//                transpose() method of corresponding matrix class.
//
// Date         : 11/19
//===========================================================================
#ifndef MATHLINEAR_H
#define MATHLINEAR_H

#include <iostream>
#include <math.h>
#include <string.h>
#include <assert.h>

namespace MATHLINEAR { // start of namespace MATHLINEAR

//===========================================================================
//
// Common math symbollic constants
//
//===========================================================================
#define ML_PI       ((float)  3.141592654f )  // Pi
#define ML_2PI      ((float)  6.283185308f )  // 2*Pi
#define ML_HLFPI    ((float)  1.570796327f )  // Pi/2 
#define ML_PIBY4    ((float)  0.785398163f )  // Pi/4
#define ML_1BYPI    ((float)  0.318309886f )  // 1/Pi
#define ML_BIGNUM   ((unsigned long) 1e6   )  // 1000000
#define ML_EPSILON  ((float)         1e-6f )  // 0.000001
#define ML_e        ((float)  2.718281828f )  // e
#define ML_LOG2e    ((float)  1.442695040f )  // log2(e)
#define ML_LOG10e   ((float)  0.434294481f )  // log10(e)
#define ML_LN2      ((float)  0.693147180f )  // ln(2)
#define ML_LN10     ((float)  2.302585092f )  // ln(10)
#define ML_SQRT2    ((float)  1.414213562f )  // Sqrt(2)
#define ML_INVSQRT2 ((float)  0.707106781f )  // 1/Sqrt(2)


//===========================================================================
//
// Some useful macros
//
//===========================================================================
#define ML_DegToRad(degree)   ( (float) (degree) * (ML_PI / 180.0f))
#define ML_RadToDeg(radian)   ( (float) (radian) * (180.0f / ML_PI))
#define ML_MIN2(x, y)         ( (x) < (y) ? (x) : (y) )
#define ML_MAX2(x, y)         ( (x) > (y) ? (x) : (y) )
#define ML_MIN3(x, y, z)      ( ML_MIN( ML_MIN((x), (y)), (z) ) )
#define ML_MAX3(x, y, z)      ( ML_MAX( ML_MAX((x), (y)), (z) ) )
#define ML_CLAMP(x, a, b)     ( (x) < (a) ? (a) : (x) > (b) ? (b) : (x) ) // clamps value 'x' in the range [a,b].
#define ML_SATURATE(x)        ( ML_CLAMP( (x), 0, 1) )					  // clamps value 'x' in the range [0,1].
#define ML_FRAND()            ( rand() / (float)RAND_MAX )                // random number within the [0,1] range.
#define ML_FRAND_RANGE(x, y)  ( (x) + ((y)-(x)) * ML_FRAND() )            // random number within [x,y] range. 
#define ML_POW2(x)            ( (a) & (a-1) == 0 ? 1 : 0     )            // tests whether the number is a power of 2      


//===========================================================================
//
// 2D Vector Class
//
//===========================================================================
struct vec2d {
	// Constructors
	inline explicit vec2d( ) { /* Empty constructor is preferred due to perf. reasons */ }
			 vec2d(float _x, float _y) { x= _x; y= _y; }
	explicit vec2d(const float *v) { assert(v); x=v[0]; y=v[1]; }
		     vec2d(const vec2d &v) { x=v.x; y=v.y; }
	
	// Unary overloaded operators
	vec2d operator ++ (int a) { return vec2d(x++, y++); }
	vec2d operator -- () { return vec2d(x--, y--); }
	vec2d operator - () const { return vec2d(-x, -y); } // vector invert

	// Binary overloaded operators
    vec2d operator + (const vec2d &v)const { return vec2d(x+v.x, y+v.y); }
    vec2d operator - (const vec2d &v)const { return vec2d(x-v.x, y-v.y); }
	// Multipl. by a scalar from rhs.
	vec2d operator * (const float scalar)const { return vec2d(scalar*x, scalar*y); } 
	// Multipl. by a scalar from lhs.
	friend   vec2d operator * (const float scalar, const vec2d &v) { return vec2d(scalar*v.x, scalar*v.y); } 
	vec2d operator / (const float scalar)const {
		 assert(scalar>ML_EPSILON);
		 float InvScalar= 1.0f / scalar;
		 return vec2d( x*InvScalar, y*InvScalar );
	}
		
	// Assignment operators
	vec2d& operator += (const vec2d &v) {
		 x += v.x;
		 y += v.y;
		 return *this;
	}
	
	vec2d& operator -= (const vec2d &v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}
	
	vec2d& operator *= (const float scalar) {
		x *= scalar;
		y *= scalar;
		return *this;
	}
	
	vec2d& operator /= (const float scalar) {
		assert(scalar);
		x /= scalar;
		y /= scalar;
		return *this;
	}
	
	vec2d& operator = (const vec2d &v) {
		x= v.x; y= v.y;
		return *this;
	}
	
	// Casting operators
	operator float * ( ) { return (float*)&x; }
	operator const float * ( )const { return (const float*)&x; }

	// Equality operators
	bool operator == (const vec2d &v) const {return ( x == v.x && y == v.y ); }
	bool operator != (const vec2d &v) const {return ( x != v.x || y != v.y ); }
	bool operator >= (const vec2d &v) const {return ( x >= v.x && y >= v.y ); }
	bool operator <= (const vec2d &v) const {return ( x <= v.x && y <= v.y ); }
	bool operator >  (const vec2d &v) const {return ( x > v.x  && y > v.y );  }
	bool operator <  (const vec2d &v) const {return ( x < v.x  && y < v.y );  } 
	
	// Set values
	void set(float vx, float vy) {x=vx; y=vy; }
	void set(const float *v){ assert(v); x=v[0]; y=v[1];}
	void setToZeroVector( ){ x = y =0.0f; }

	// console I/O 
	inline friend std::ostream& operator << (std::ostream &out, const vec2d &v) {
		return out << "(" << v.x << "," << v.y << ")" <<'\n';
	}
	
	inline friend std::istream& operator >> (std::istream &in, vec2d &v) {
		return in >> v.x >> v.y;
	}
	
	// Dot product
	inline float dot(const vec2d &v) const { return x*v.x + y*v.y; } 

	// Cross product
	inline float cross(const vec2d &v) const { return x*v.y - y*v.x; } 
	
	// Magnitude
	inline float length() const { return sqrtf( dot(*this) ); } 

	// Magnitude squared
	inline float lengthSqr()const { return dot(*this); }

	// Distance between two vectors
	inline float distance(const vec2d &v)const {
		return ((*this) - v).length();
	}

	// Normalize
	inline void normalize() {
		float mag = this->length();
		if(mag>ML_EPSILON) {
			this->x *= 1.0f/mag;
			this->y *= 1.0f/mag;
		}
	}

	// Return angle(radian) between two vectors
	inline float angleBtw(const vec2d &v) {
		float invDen = this->length() *  v.length();
		assert(invDen>ML_EPSILON);
		// clamp the arccos range to [-1,1]
		return acos(ML_CLAMP( dot(v)/invDen, -1.0f, 1.0f));
	}

	// Reflection vector R = v2 - 2*V1*dot(v1*v2)
	inline vec2d reflection(const vec2d &v) const {
		return vec2d( (*this)- 2.0f*v*dot(v) );
	}

	// Linear interpolation between two vectors i.e. V1*(1-s) + V2*s
	inline vec2d Lerp(const vec2d &v, const float s) const {
		return vec2d( v*s + (1 - s)*(*this) );
	}
	
	// 2D Coordinates
	union {
		struct {
			float x;
			float y;
		};
		float v[2];
	};

};


//===========================================================================
//
// 3D Vector
//
//===========================================================================
struct vec3d {
	// Constructors
	explicit vec3d( ) { /* Empty constructor is preferred due to perf. reasons */ }
			 vec3d(float _x, float _y, float _z) { x= _x; y= _y; z=_z; }
			 vec3d(const vec2d &v, float scalar) : x(v.x), y(v.y), z(scalar) { }
	explicit vec3d(const float *v) { assert(v); x=v[0]; y=v[1]; z=v[2]; }
		     vec3d(const vec3d &v) { x=v.x; y=v.y; z=v.z; }
	
	// Unary overloaded operators
	vec3d operator ++ (int a) { return vec3d(x++, y++, z++); }
	vec3d operator -- () { return vec3d(x--, y--, z--); }
	vec3d operator - () const { return vec3d(-x, -y, -z); } // vector invert

	// Binary overloaded operators
    vec3d operator + (const vec3d &v)const { return vec3d(x+v.x, y+v.y, z+v.z); }
    vec3d operator - (const vec3d &v)const { return vec3d(x-v.x, y-v.y, z-v.z); }
	// Multipl. by a scalar from rhs.
	vec3d operator * (const float scalar)const { return vec3d(scalar*x, scalar*y, scalar*z); } 
	// Multipl. by a scalar from lhs.
	friend   vec3d operator * (const float scalar, const vec3d &v) { return vec3d(scalar*v.x, scalar*v.y, scalar*v.z); } 
	vec3d operator / (const float scalar)const {
		 assert(scalar);
		 float InvScalar = 1.0f / scalar;
		 return vec3d(x*InvScalar, y*InvScalar, z*InvScalar);
	}
		
	// Assignment operators
	vec3d& operator += (const vec3d &v) {
		 x += v.x;
		 y += v.y;
		 z += v.z;
		 return *this;
	}
	
	vec3d& operator -= (const vec3d &v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
	
	vec3d& operator *= (const float scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}
	
	vec3d& operator /= (const float scalar) {
		assert(scalar);
		x /= scalar;
		y /= scalar;
		z /= scalar;
		return *this;
	}
	
	vec3d& operator = (const vec3d &v) {
		x= v.x; y= v.y; z= v.z;
		return *this;
	}
	
	// Casting operators
	operator float * () { return (float*)&x; }
	operator const float * ()const { return (const float*)&x; }

	// Equality operators
	bool operator == (const vec3d &v) const {return ( x == v.x && y == v.y && z==v.z ); }
	bool operator != (const vec3d &v) const {return ( x != v.x || y != v.y || z != v.z); }
	bool operator >= (const vec3d &v) const {return ( x >= v.x && y >= v.y && z >= v.z); }
	bool operator <= (const vec3d &v) const {return ( x <= v.x && y <= v.y && y <= v.z); }
	bool operator >  (const vec3d &v) const {return ( x > v.x && y > v.y && z > v.z);   }
	bool operator <  (const vec3d &v) const {return ( x < v.x && y < v.y && z <v.z );   } 
	
	// Set values
	void set( float vx, float vy, float vz) {x=vx; y=vy; z=vz; }
	void set(const float *v){ assert(v); x=v[0]; y=v[1]; z=v[2];}
	void setToZeroVector( ){ x = y = z = 0.0f; }

	// console I/O 
	inline friend std::ostream& operator << (std::ostream &out, const vec3d &v) {
		return out << "(" << v.x << "," << v.y << "," << v.z << ")" << '\n';
	}
	
	inline friend std::istream& operator >> (std::istream &in, vec3d &v) {
		return in >> v.x >> v.y >> v.z;
	}
	
	// Dot product
	inline float dot(const vec3d &v) const { return x*v.x + y*v.y + z*v.z; } 

	// Cross product
	inline vec3d cross(const vec3d &v) const { 
		return vec3d(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
	}

	// Magnitude
	inline float length() const { return sqrtf( dot(*this) ); } 

	// Magnitude squared
	inline float lengthSqr()const { return dot(*this); }

	// Distance between two vectors
	inline float distance(const vec3d &v)const {
		return ((*this) - v).length();
	}

	// Normalize
	inline void normalize() {
		float mag = this->length();
		if(mag>ML_EPSILON) {
			this->x *= 1.0f/mag;
			this->y *= 1.0f/mag;
			this->z *= 1.0f/mag;

		}
	}

	// Return angle(radian) between two vectors
	inline float angleBtw(const vec3d &v) {
		float invDen = this->length() *  v.length();
		assert(invDen>ML_EPSILON);
		// clamp the arccos range to [-1,1]
		return acos( ML_CLAMP( dot(v)/invDen, -1.0f, 1.0f) );
	}

	// Reflection vector R = v2 - 2*v1*dot(v1*v2)
	inline vec3d reflection(const vec3d &v) const {
		return vec3d( (*this)- 2.0f*v*dot(v) );
	}

	// Linear interpolation between vectors
	inline vec3d Lerp(const vec3d &v, const float s) const {
		return vec3d( (*this) + (v - (*this))*s);
	}
	
	// 3D Coordinates
	union {
		struct {
			float x;
			float y;
			float z;
		};
		float v[3];
	};

};


//===========================================================================
//
// 4D Vector
//
//===========================================================================
struct vec4d {
	// Constructors
	explicit vec4d( ) { /* Empty constructor is preferred due to perf. reasons */ }
		     vec4d(float _x, float _y, float _z, float _w=1.0f) : x(_x), y(_y), z(_z), w(_w) {  }
			 vec4d(const vec2d &v1, const vec2d &v2) : x(v1.x), y(v1.y), z(v2.x), w(v2.y) { }
			 vec4d(const vec3d &v, float _w) : x(v.x), y(v.y), z(v.z), w(_w) { } 
	explicit vec4d(const float *v) { assert(v); x=v[0]; y=v[1]; z=v[2]; w=v[3]; }
		     vec4d(const vec4d &v) { x=v.x; y=v.y; z=v.z; w=v.w; }
	
	// Unary overloaded operators
	vec4d operator ++ (int a) { return vec4d(x++, y++, z++, w++); }
	vec4d operator -- () { return vec4d(x--, y--, z--, w--); }
	vec4d operator - () const { return vec4d(-x, -y, -z, -w); } // vector invert

	// Binary overloaded operators
    vec4d operator + (const vec4d &v)const { return vec4d(x+v.x, y+v.y, z+v.z, w+v.w); }
    vec4d operator - (const vec4d &v)const { return vec4d(x-v.x, y-v.y, z-v.z, w-v.w ); }
	// Multipl. by a scalar from rhs.
	vec4d operator * (const float scalar)const { return vec4d(scalar*x, scalar*y, scalar*z, scalar*w); } 
	// Multipl. by a scalar from lhs.
	friend   vec4d operator * (const float scalar, const vec4d &v) { return vec4d(scalar*v.x, scalar*v.y, scalar*v.z, scalar*v.w); } 
	vec4d operator / (const float scalar)const {
		 assert(scalar>ML_EPSILON);
		 float InvScalar = 1.0f / scalar;
		 return vec4d(x*InvScalar, y*InvScalar, z*InvScalar, w*InvScalar);
	}
		
	// Assignment operators
	vec4d& operator += (const vec4d &v) {
		 x += v.x;
		 y += v.y;
		 z += v.z;
		 w += v.w;
		 return *this;
	}
	
	vec4d& operator -= (const vec4d &v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		w -= v.w;
		return *this;
	}
	
	vec4d& operator *= (const float scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
		return *this;
	}
	
	vec4d& operator /= (const float scalar) {
		assert(scalar);
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;
		return *this;
	}
	
	vec4d& operator = (const vec4d &v) {
		x= v.x; y= v.y; z= v.z; w = v.w;
		return *this;
	}
	
	// Casting operators
	operator float * () { return (float*)&x; }
	operator const float * ()const { return (const float*)&x; }

	// Equality operators
	bool operator == (const vec4d &v) const {return ( x == v.x && y == v.y && z==v.z && w == v.w );  }
	bool operator != (const vec4d &v) const {return ( x != v.x || y != v.y || z != v.z || w != v.w); }
	bool operator >= (const vec4d &v) const {return ( x >= v.x && y >= v.y && z >= v.z && w >= v.w); }
	bool operator <= (const vec4d &v) const {return ( x <= v.x && y <= v.y && y <= v.z && w <= v.w); }
	bool operator >  (const vec4d &v) const {return ( x > v.x  && y > v.y  && z > v.z  && w > v.w);  }
	bool operator <  (const vec4d &v) const {return ( x < v.x  && y < v.y  && z < v.z  && w < v.w);  } 
	
	// Set values
	void set( float vx, float vy, float vz, float vw) {x=vx; y=vy; z=vz; w=vw; }
	void set(const float *v){ assert(v); x=v[0]; y=v[1]; z=v[2]; w=v[3];}
	void setToZeroVector( ){ x= y= z= w= 0.0f; }

	// console I/O 
	inline friend std::ostream& operator << (std::ostream &out, const vec4d &v) {
		return out << "(" << v.x << "," << v.y << "," << v.z << "," << v.w << ")" <<'\n';
	}
	
	inline friend std::istream& operator >> (std::istream &in, vec4d &v) {
		return in >> v.x >> v.y >> v.z >> v.w;
	}
	
	// Dot product
	inline float dot(const vec4d &v) const { return x*v.x + y*v.y + z*v.z + w*v.w; } 

	// Cross product
	//inline vec4d cross(const vec4d &v) const { return vec4d(); }

	// Magnitude
	inline float length() const { return sqrtf( dot(*this) ); } 

	// Magnitude squared
	inline float lengthSqr()const { return dot(*this); }

	// Distance between two vectors
	inline float distance(const vec4d &v)const {
		return ((*this) - v).length();
	}

	// Normalize
	inline void normalize() {
		float mag = this->length();
		if(mag>ML_EPSILON) {
			this->x *= 1.0f/mag;
			this->y *= 1.0f/mag;
			this->z *= 1.0f/mag;
			this->w *= 1.0f/mag;
		}
	}

	// Return angle(radian) between two vectors
	inline float angleBtw(const vec4d &v) {
		float invDen = this->length() *  v.length();
		assert(invDen>ML_EPSILON);
		// clamp the arccos range to [-1,1]
		return acos( ML_CLAMP( dot(v)/invDen, -1.0f, 1.0f) );
	}

	// Reflection vector R = v2 - 2*V1*dot(v1*v2)
	inline vec4d reflection(const vec4d &v) const {
		return vec4d( (*this)- 2.0f*v*dot(v) );
	}

	// Linear interpolation between vectors
	inline vec4d Lerp(const vec4d &v, const float s) const {
		return vec4d(*this + (v - *this)*s);
	}
	
	// 4D Coordinates
	union {
		struct {
			float x;
			float y;
			float z;
			float w;
		};
		float v[4];
	};

};

//===========================================================================
//
// 3x3 Matrix Class
//
//===========================================================================
struct mat3x3 {
	// Constructors
	explicit mat3x3() { /* Empty constructor is preferred due to perf. reasons */ }
			 mat3x3( float _m11, float _m12, float _m13,
					 float _m21, float _m22, float _m23,
					 float _m31, float _m32, float _m33 ) {
			         
				m11= _m11; m12= _m12; m13= _m13;
				m21= _m21; m22= _m22; m23= _m23;
				m31= _m31; m32= _m32; m33= _m33;
			}
			 
	explicit mat3x3(const float *marray) {
				 assert(marray);
				 memcpy(this, marray, sizeof(float)*9);
	}

			 mat3x3(const vec3d &r1, const vec3d &r2, const vec3d &r3) {
				 
				 memcpy(&m11, &r1, sizeof(float)*3);
				 memcpy(&m21, &r2, sizeof(float)*3);
				 memcpy(&m31, &r3, sizeof(float)*3);
			 }
			 mat3x3(const mat3x3 &mat) {
				 
				m11= mat.m11; m12= mat.m12; m13= mat.m13;
				m21= mat.m21; m22= mat.m22; m23= mat.m23;
				m31= mat.m31; m32= mat.m32; m33= mat.m33;
			 }
				 
	
	// Set Matrix row values
	inline void set( float _m11, float _m12, float _m13,
					 float _m21, float _m22, float _m23,
					 float _m31, float _m32, float _m33 ) {
				
		m11= _m11; m12= _m12; m13= _m13;
		m21= _m21; m22= _m22; m23= _m23;
		m31= _m31; m32= _m32; m33= _m33;
	}
	inline void set(const vec3d &r1, const vec3d &r2, const vec3d &r3){
		memcpy(&m11, &r1, sizeof(vec3d));
		memcpy(&m21, &r2, sizeof(vec3d));
		memcpy(&m31, &r3, sizeof(vec3d));
	}
	inline void setFromArray(const float *m) {
		memcpy(this, m, sizeof(float)*9);
	}
	inline vec3d getRow(unsigned row) const {
		assert(row>=1 && row<=3);
		return vec3d( m[3*(row-1)+0], m[3*(row-1)+1], m[3*(row-1)+2] );
	}
	inline vec3d getCol(unsigned col) const {
		assert(col>=1 && col<=3);
		return vec3d( m[3*0+col-1], m[3*1+col-1], m[3*2+col-1] );
	}
	inline void identity() {
		
		m11 = m22 = m33 =  1.0f;
		m12 = m13 = m21 = m23 = m31 = m32 = 0.0f;
	}
	
	// Unary minus operator
	inline mat3x3 operator - () const {

		return mat3x3(-m11, -m12, -m13, 
			          -m21, -m22, -m23,
					  -m31, -m32, -m33 );
	}
	
	// Binary Operators
	inline mat3x3 operator + (const mat3x3 &mat) const {

		return mat3x3( m11+mat.m11, m12+mat.m12, m13+mat.m13, 
					   m21+mat.m21, m22+mat.m22, m23+mat.m23, 
					   m31+mat.m31, m32+mat.m32, m33+mat.m33 );
	}
	
	inline mat3x3 operator - (const mat3x3 &mat) const {
		
		return mat3x3( m11-mat.m11, m12-mat.m12, m13-mat.m13, 
					   m21-mat.m21, m22-mat.m22, m23-mat.m23, 
					   m31-mat.m31, m32-mat.m32, m33-mat.m33 );
	}
	
	// Multiplication with a scalar from rhs. 
	inline mat3x3 operator * (const float scalar) const {
		
		return mat3x3( m11*scalar, m12*scalar, m13*scalar, 
					   m21*scalar, m22*scalar, m23*scalar, 
					   m31*scalar, m32*scalar, m33*scalar );
	
	}
	
	// Multiplication with a scalar from lhs.
	friend mat3x3 operator * (const float scalar, const mat3x3 &mat) {

		return mat3x3( scalar*mat.m11, scalar*mat.m12, scalar*mat.m13, 
					   scalar*mat.m21, scalar*mat.m22, scalar*mat.m23, 
					   scalar*mat.m31, scalar*mat.m32, scalar*mat.m33 );
	}
		
	// matrix multiplication
	mat3x3 operator * (const mat3x3 &mat) const {

		mat3x3 temp;
		for(int i=0;i<3;++i) {
			for(int j=0;j<3;++j) {
				temp.m[3*i+j] = 0.0f;
				for(int k=0;k<3;++k)
					temp.m[i*3+j] += m[3*i+k] * mat.m[3*k+j];
			}
		}
		return temp;
	}

	// Matrix-Vector multiplication: M*v
	vec3d operator * (const vec3d &v) {

		return vec3d( m11*v.x + m12*v.y + m13*v.z,
			          m21*v.x + m22*v.y + m23*v.z,
					  m31*v.x + m32*v.y + m33*v.z );
	}
	
	// Vector-Matrix multiplication: v*M
	inline friend vec3d operator * (const vec3d &v, const mat3x3 &m) {
		
		return vec3d( v.x*m.m11 + v.y*m.m21 + v.z*m.m31,
			          v.x*m.m12 + v.y*m.m22 + v.z*m.m32,
					  v.x*m.m13 + v.y*m.m23 + v.z*m.m33 );
	}

	mat3x3 operator / (const float scalar) const {
		
		assert(scalar);
		float invs = 1.0f/scalar;
		return mat3x3( m11*invs, m12*invs, m13*invs, 
					   m21*invs, m22*invs, m23*invs, 
					   m31*invs, m32*invs, m33*invs );
		
	}
		
	// Assignment operators
	mat3x3& operator += (const mat3x3 &mat) {

		return *this = *this + mat;
	}
	
	mat3x3& operator -= (const mat3x3 &mat) {
		
		return *this = *this - mat;
	}
	
	mat3x3& operator *= (const float scalar) {
		
		return *this = *this * scalar;
	}

	mat3x3& operator *= (const mat3x3 &mat) {

		return *this = *this * mat;
	}
	
	mat3x3& operator /= (const float scalar) {
		
		assert(scalar);
		return *this = *this / scalar;
	}
		
	// Casting operators
	operator       float* ( )       { return (float*)&m11;       }
	operator const float* ( ) const { return (const float*)&m11; }

	// array index operators
	float& operator ( ) ( unsigned int Row, unsigned int Col ) {

		// Check bounds first
		assert( (Row <= 3  && Row >= 1) && (Col <= 3 && Col >= 1) ); 
		return *( &this->m11 + (Row-1)*3 + (Col-1) );
	}
	
	float  operator ( ) ( unsigned int Row, unsigned int Col ) const {
		// Check bounds first
		assert( (Row <= 3  && Row >= 1) && (Col <= 3 && Col >= 1) ); 
		return *( &this->m11 + (Row-1)*3 + (Col-1) );
	}

	// array subscript operator
	float operator [] (unsigned int i) {
		assert( i>=1 && i<=9 );
		return *(&this->m11 + i-1);
	}
	
	// Equality operators
	bool operator == (const mat3x3 &mat) const {
		return  0 == memcmp( this, &mat, sizeof(mat3x3));
	}
	
	bool operator != (const mat3x3 &mat) const {
		return 0 != memcmp( this, &mat, sizeof(mat3x3));
	}
	
	// Prints matrix data to the default output stream
	friend std::ostream& operator << ( std::ostream &out, const mat3x3 &mat ) {
		return out <<
		"|" << mat.m11 <<"  "<< mat.m12<<"  "<<mat.m13<<"|\n"
		"|" << mat.m21 <<"  "<< mat.m22<<"  "<<mat.m23<<"|\n"
		"|" << mat.m31 <<"  "<< mat.m32<<"  "<<mat.m33<<"|\n";
	}
	
	float trace( ) const { return m11 + m22 + m33; }
	
	void transposeSelf( )  {
		float temp;
		temp = m21; m21=m12; m12=temp;
		temp = m31; m31=m13; m13=temp;
		temp = m32; m32=m23; m23=temp;
		
	}
	
	mat3x3 getTranspose() const {
		return mat3x3(m11, m21, m31,
			          m12, m22, m32,
					  m13, m23, m33 ); 
	}

	void inverseSelf( ) {

		float det = determinant();
		assert( fabs(det) > ML_EPSILON && "Error: Singular matrix...");
		float invDet =  1.0f / det;

		// Calculate cofactor terms
		float cof11 = m22*m33 - m23*m32;
		float cof21 = m23*m31 - m21*m33;
		float cof31 = m21*m32 - m22*m31;
		float cof12 = m13*m32 - m12*m33;
		float cof22 = m11*m33 - m13*m31;
		float cof32 = m12*m31 - m11*m32;
		float cof13 = m12*m23 - m13*m22;
		float cof23 = m13*m21 - m11*m23;
		float cof33 = m11*m22 - m12*m21;

		m11= cof11*invDet; m12= cof12*invDet; m13= cof13*invDet;
		m21= cof21*invDet; m22= cof22*invDet; m23= cof23*invDet;
		m31= cof31*invDet; m32= cof32*invDet; m33= cof33*invDet;

	}
	
	// returns an inverse matrix copy w/o changing original matrix data.
	mat3x3 getInverse( ) const {
		float det = determinant();
		assert( fabs(det) > ML_EPSILON && "Error: Singular matrix...");
		float invDet =  1.0f / det;

		// Calculate cofactor terms
		float cof11 = m22*m33 - m23*m32;
		float cof21 = m23*m31 - m21*m33;
		float cof31 = m21*m32 - m22*m31;
		float cof12 = m13*m32 - m12*m33;
		float cof22 = m11*m33 - m13*m31;
		float cof32 = m12*m31 - m11*m32;
		float cof13 = m12*m23 - m13*m22;
		float cof23 = m13*m21 - m11*m23;
		float cof33 = m11*m22 - m12*m21;
		
		return mat3x3( cof11*invDet, cof12*invDet, cof13*invDet,
			           cof21*invDet, cof22*invDet, cof23*invDet,
					   cof31*invDet, cof32*invDet, cof33*invDet );
		
	}
	
	float determinant( ) const {
		
		// calculate sub determinants
		float det11 =    m22*m33 - m23*m32;
		float det12 =  -(m21*m33 - m23*m31);
		float det13 =    m21*m32 - m22*m31;

		return m11*det11 + m12*det12 + m13*det13;
	}
	
	// Transformations
	void makeTranslateMatrix(float x, float y, float z) {
		identity();
		m13 = x; m23 = y; m33 = z;
	}
	
	void makeTranslateMatrix(const vec3d &t) {
		makeTranslateMatrix(t.x, t.y, t.z);
	}

	void makeRotXMatrix(float angle) {
		float c = cos(angle), s = sin(angle);
		identity();
		m22 = c; m23 = -s;
		m32 = s; m33 = c;
	}
	
	void makeRotYMatrix(float angle) {
		float c = cos(angle), s = sin(angle);
		identity();
		m11 = c;  m13 = s;
		m31 = -s; m33 = c;
	}
	
	void makeRotZMatrix(float angle) {
		float c = cos(angle), s = sin(angle);
		identity();
		m11 = c; m12 = -s;
		m21 = s; m22 = c;
	}
	
	void makeRotAxis(float angle, const vec3d &v) {
		
		float s=sin(ML_DegToRad(angle)), c=cos(ML_DegToRad(angle)), onemc=1.0f-c;
		float xy, xz, yz;
		vec3d axis(v);

		// guarantee unit rot. axis
		if(1.0f != axis.length()) axis.normalize();

		xy = axis.x*axis.y;
		xz = axis.x*axis.z;
		yz = axis.y*axis.z;
		m11 = axis.x*axis.x*onemc + c;
		m12 = xy*onemc - axis.z*s;
		m13 = xz*onemc + axis.y*s;
		m21 = xy*onemc + axis.z*s;
		m22 = axis.y*axis.y*onemc + c;
		m23 = yz*onemc - axis.x*s;
		m31 = xz*onemc - axis.y*s;
		m32 = yz*onemc + axis.x*s;
		m33 = axis.z*axis.z*onemc + c;
	}
	
	void makeScaleMatrix( float sx, float sy, float sz) {
		identity();
		m11 = sx; m22 = sy; m33 = sz;
	}

	// Query methods
	bool isSymmetric( ) const { return *this == getTranspose(); }
	
	bool isDiagonal ( ) const {
		return (  m12== 0.0f && m13== 0.0f && m21== 0.0f &&
			      m23== 0.0f && m31== 0.0f && m32== 0.0f );
	}
	
	bool isIdentity ( ) const {
		return ( m11==1.0f && m22==1.0f && m33==1.0f ) &&
	           ( m12==0.0f && m13==0.0f && m21==0.0f ) &&
			   ( m23==0.0f && m31==0.0f && m32==0.0f); 
	}
	
	union {
		struct {
			float m11, m12, m13;
			float m21, m22, m23;
			float m31, m32, m33;
		};
		float m[9];
	};
};


//===========================================================================
//
// 4x4 Matrix Class
//
//===========================================================================
struct mat4x4 {
	// Constructors
	explicit mat4x4( ) { /* Empty constructor is preferred due to perf. reasons */ }
			 mat4x4(float _m11, float _m12, float _m13, float _m14,
					float _m21, float _m22, float _m23, float _m24,
				    float _m31, float _m32, float _m33, float _m34,
					float _m41=0.0f, float _m42=0.0f, float _m43=0.0f, float _m44=1.0f ) {

					m11= _m11; m12= _m12; m13= _m13; m14= _m14;
					m21= _m21; m22= _m22; m23= _m23; m24= _m24;
					m31= _m31; m32= _m32; m33= _m33; m34= _m34;
					m41= _m41; m42= _m42; m43= _m43; m44= _m34;
			 }
			 
			 mat4x4(const vec4d &r1, const vec4d &r2, 
				    const vec4d &r3, const vec4d &r4 ) {

					memcpy(&m11, &r1, sizeof(float)*4);
					memcpy(&m21, &r2, sizeof(float)*4);
					memcpy(&m31, &r3, sizeof(float)*4);
					memcpy(&m41, &r4, sizeof(float)*4);
			 }

	explicit mat4x4(const float *m) {
			assert(m);
			memcpy(&m11, m, sizeof(float)*16);
	}
			mat4x4(const mat4x4 &mat) {

				memcpy(this, &mat, sizeof(float)*16);
			}

	
	// Set Matrix row values
	inline void set( float _m11, float _m12, float _m13, float _m14,
				     float _m21, float _m22, float _m23, float _m24,
					 float _m31, float _m32, float _m33, float _m34,
					 float _m41=0.0f, float _m42=0.0f, float _m43=0.0f, float _m44=1.0f ) {

		m11= _m11; m12= _m12; m13= _m13; m14= _m14;
		m21= _m21; m22= _m22; m23= _m23; m24= _m24;
		m31= _m31; m32= _m32; m33= _m33; m34= _m34;
		m41= _m41; m42= _m42; m43= _m43; m44= _m44;
	}
	inline void set ( const vec4d &r1, const vec4d &r2,
		              const vec4d &r3, const vec4d &r4 ) {
		memcpy(&m11, &r1, sizeof(vec4d));
		memcpy(&m21, &r2, sizeof(vec4d));
		memcpy(&m31, &r3, sizeof(vec4d));
		memcpy(&m41, &r4, sizeof(vec4d));
	}

	inline void set(const float *m) {
		assert(m);
		memcpy(this, m, sizeof(float)*16);
	}

	inline vec4d getRow(unsigned row) const {
		assert(row>=1 && row<=1);
		return vec4d(m[4*(row-1)+0],m[4*(row-1)+1],m[4*(row-1)+2],m[4*(row-1)+3]);
	}

	inline vec4d getCol(unsigned col) const {
		assert(col>=1 && col<=1);
		return vec4d(m[4*0+col-1],m[4*1+col-1],m[4*2+col-1],m[4*3+col-1]);
	}

	inline void identity() {
		memset(this, 0, sizeof(float)*16);
		m11= m22= m33= m44= 1.0f;
    }
	
	// Unary minus operator
	inline mat4x4 operator - () const {

		return mat4x4(-m11, -m12, -m13, -m14 
			          -m21, -m22, -m23, -m24,
					  -m31, -m32, -m33, -m34,
					  -m41, -m42, -m43, -m44 );
	}
	
	// Binary Operators
	inline mat4x4 operator + (const mat4x4 &mat) const {

		return mat4x4( m11+mat.m11, m12+mat.m12, m13+mat.m13,m14+mat.m14, 
					   m21+mat.m21, m22+mat.m22, m23+mat.m23,m24+mat.m24, 
					   m31+mat.m31, m32+mat.m32, m33+mat.m33,m34+mat.m34,
					   m41+mat.m41, m42+mat.m42, m43+mat.m43,m44+mat.m44 );
	}
	
	inline mat4x4 operator - (const mat4x4 &mat) const {
		
		return mat4x4( m11-mat.m11, m12-mat.m12, m13-mat.m13,m14-mat.m14, 
					   m21-mat.m21, m22-mat.m22, m23-mat.m23,m24-mat.m24, 
					   m31-mat.m31, m32-mat.m32, m33-mat.m33,m34-mat.m34,
					   m41-mat.m41, m42-mat.m42, m43-mat.m43,m44-mat.m44 );
	}
	
	// Multiplication with a scalar from rhs. 
	inline mat4x4 operator * (const float scalar) const {
		
		return mat4x4( m11*scalar, m12*scalar, m13*scalar, m14*scalar, 
					   m21*scalar, m22*scalar, m23*scalar, m24*scalar,
					   m31*scalar, m32*scalar, m33*scalar, m34*scalar,
					   m41*scalar, m42*scalar, m43*scalar, m44*scalar );
	
	}
	
	// Multiplication with a scalar from lhs.
	friend mat4x4 operator * (const float scalar, const mat4x4 &mat) {

		return mat4x4( scalar*mat.m11,scalar*mat.m12,scalar*mat.m13,scalar*mat.m14, 
					   scalar*mat.m21,scalar*mat.m22,scalar*mat.m23,scalar*mat.m24,
					   scalar*mat.m31,scalar*mat.m32,scalar*mat.m33,scalar*mat.m34,
					   scalar*mat.m31,scalar*mat.m32,scalar*mat.m33,scalar*mat.m44 );
					   
	}
		
	// matrix multiplication
	mat4x4 operator * (const mat4x4 &mat) const {

		mat4x4 temp;
		for(int i=0;i<4;++i) {
			for(int j=0;j<4;++j) {
				temp.m[4*i+j] = 0.0f;
				for(int k=0;k<4;++k)
					temp.m[i*4+j] += m[4*i+k] * mat.m[4*k+j];
			}
		}
		return temp;
	}

	// Matrix-Vector multiplication: M*v
	vec4d operator * (const vec4d &v) {

		return vec4d( m11*v.x + m12*v.y + m13*v.z + m14*v.w,
			          m21*v.x + m22*v.y + m23*v.z + m24*v.w, 
					  m31*v.x + m32*v.y + m33*v.z + m34*v.w,
					  m41*v.x + m42*v.y + m43*v.z + m44*v.w );
	}
	
	// Vector-Matrix multiplication: v*M
	inline friend vec4d operator * (const vec4d &v, const mat4x4 &m) {
		
		return vec4d( v.x*m.m11 + v.y*m.m21 + v.z*m.m31 + v.w*m.m41,
			          v.x*m.m12 + v.y*m.m22 + v.z*m.m32 + v.w*m.m42,
					  v.x*m.m13 + v.y*m.m23 + v.z*m.m33 + v.w*m.m43,
					  v.x*m.m14 + v.y*m.m24 + v.z*m.m34 + v.w*m.m44 );
	}

	mat4x4 operator / (const float scalar) const {
		
		assert(scalar);
		float invs = 1.0f/scalar;
		return mat4x4( m11*invs, m12*invs, m13*invs, m14*invs,  
					   m21*invs, m22*invs, m23*invs, m24*invs,
					   m31*invs, m32*invs, m33*invs, m34*invs,
					   m41*invs, m42*invs, m43*invs, m44*invs );
		
	}
		
	// Assignment operators
	mat4x4& operator += (const mat4x4 &mat) {

		return *this = *this + mat;
	}
	
	mat4x4& operator -= (const mat4x4 &mat) {
		
		return *this = *this - mat;
	}
	
	mat4x4& operator *= (const float scalar) {
		
		return *this = *this * scalar;
	}

	mat4x4& operator *= (const mat4x4 &mat) {

		return *this = *this * mat;
	}
	
	mat4x4& operator /= (const float scalar) {
		
		assert(scalar);
		return *this = *this / scalar;
	}
		
	// Casting operators
	operator       float* ( )       { return (float*)&m11;       }
	operator const float* ( ) const { return (const float*)&m11; }

	// array index operators
	float& operator ( ) ( unsigned Row, unsigned Col ) {

		// Check bounds first
		assert( (Row <= 4  && Row >= 1) && (Col <= 4 && Col >= 1) ); 
		return *( &this->m11 + (Row-1)*4 + (Col-1) );
	}
	
	float  operator ( ) ( unsigned Row, unsigned Col ) const {
		// Check bounds first
		assert( (Row <= 4  && Row >= 1) && (Col <= 4 && Col >= 1) ); 
		return *( &this->m11 + (Row-1)*4 + (Col-1) );
	}

	// array subscript operator
	float operator [] (unsigned i) {
		assert( i>=1 && i<=16 );
		return *(&this->m11 + i-1);
	}
	
	// Equality operators
	bool operator == (const mat4x4 &mat) const {
		return  0 == memcmp( this, &mat, sizeof(mat4x4));
	}
	
	bool operator != (const mat4x4 &mat) const {
		return 0 != memcmp( this, &mat, sizeof(mat4x4));
	}
	
	// Prints matrix data to the default output stream
	friend std::ostream& operator << ( std::ostream &out, const mat4x4 &mat ) {
		return out <<
			"|" << mat.m11 <<"  "<< mat.m12<<"  "<<mat.m13<< " "<< mat.m14 << "|\n"
			"|" << mat.m21 <<"  "<< mat.m22<<"  "<<mat.m23<< " "<< mat.m24 << "|\n"
			"|" << mat.m31 <<"  "<< mat.m32<<"  "<<mat.m33<<" "<<  mat.m34 << "|\n"
			"|" << mat.m41 <<"  "<< mat.m42<<"  "<<mat.m43<<" "<<  mat.m44 << "|\n";
	}
	
	float trace( ) const { return m11 + m22 + m33 + m44; }
	
	void  transposeSelf( )  {
		float temp;
		temp = m21; m21=m12; m12=temp;
		temp = m31; m31=m13; m13=temp;
		temp = m32; m32=m23; m23=temp;
		temp = m41; m41=m14; m14=temp;
		temp = m42; m42=m24; m24=temp;
		temp = m43; m43=m34; m34=temp;
	
	}

	mat4x4 getTranspose()const {
		return mat4x4( m11, m21, m31, m41,
			           m12, m22, m32, m42,
					   m13, m23, m33, m43,
					   m14, m24, m34, m44 );
	}

	void inverseSelf( ) {

		float det = determinant(), invDet;
		assert(fabs(det)>ML_EPSILON && "Error: Singular matrix");
		invDet =  1.0f / det;

		// Calculate cofactor terms
		float cof11= m22*m33*m44+m23*m34*m42+m24*m32*m43-m22*m34*m43-m23*m32*m44-m24*m33*m42;
		float cof12= m12*m34*m43+m13*m32*m44+m14*m33*m42-m12*m33*m44-m13*m34*m42-m14*m32*m43; 
		float cof13= m12*m23*m44+m13*m24*m42+m14*m22*m43-m12*m24*m43-m13*m22*m44-m14*m23*m42; 
		float cof14= m12*m24*m33+m13*m22*m34+m14*m23*m32-m12*m23*m34-m13*m24*m32-m14*m22*m33;
		
		float cof21= m21*m34*m43+m23*m31*m44+m24*m33*m41-m21*m33*m44-m23*m34*m41-m24*m31*m43;
		float cof22= m11*m33*m44+m13*m34*m41+m14*m31*m43-m11*m34*m43-m13*m31*m44-m14*m33*m41;
		float cof23= m11*m24*m43+m13*m21*m44+m14*m23*m41-m11*m23*m44-m13*m24*m41-m14*m21*m43;
		float cof24= m11*m23*m34+m13*m24*m31+m14*m21*m33-m11*m24*m33-m13*m21*m34-m14*m23*m31;
		
		float cof31= m21*m32*m44+m22*m34*m41+m24*m31*m42-m21*m34*m42-m22*m31*m44-m24*m32*m41;
		float cof32= m11*m34*m42+m12*m31*m44+m14*m32*m41-m11*m32*m44-m12*m34*m41-m14*m31*m42;
		float cof33= m11*m22*m44+m12*m24*m41+m14*m21*m42-m11*m24*m42-m12*m21*m44-m14*m22*m41;
		float cof34= m11*m24*m32+m12*m21*m34+m14*m22*m31-m11*m22*m34-m12*m24*m31-m14*m21*m32;
		
		float cof41= m21*m33*m42+m22*m31*m43+m23*m32*m41-m21*m32*m43-m22*m33*m41-m23*m31*m42;
		float cof42= m11*m32*m43+m12*m33*m41+m13*m31*m42-m11*m33*m42-m12*m31*m43-m13*m32*m41;
		float cof43= m11*m23*m42+m12*m21*m43+m13*m22*m41-m11*m22*m43-m12*m23*m41-m13*m21*m42;
		float cof44= m11*m22*m33+m12*m23*m31+m13*m21*m32-m11*m23*m32-m12*m21*m33-m13*m22*m31;
		
		m11= cof11*invDet; m12= cof12*invDet; m13= cof13*invDet; m14= cof14*invDet;
		m21= cof21*invDet; m22= cof22*invDet; m23= cof23*invDet; m24= cof24*invDet;
		m31= cof31*invDet; m32= cof32*invDet; m33= cof33*invDet; m34= cof34*invDet;
		m41= cof41*invDet; m42= cof42*invDet; m43= cof43*invDet; m44= cof44*invDet;
		 
	}

	// returns an inverse matrix copy w/o changing the original matrix data.
	mat4x4 getInverse() const {
		
		float det = determinant(), invDet;
		assert(fabs(det)>ML_EPSILON && "Error: Singular matrix");
		invDet =  1.0f / det;

		// Calculate cofactor terms
		float cof11= m22*m33*m44+m23*m34*m42+m24*m32*m43-m22*m34*m43-m23*m32*m44-m24*m33*m42;
		float cof12= m12*m34*m43+m13*m32*m44+m14*m33*m42-m12*m33*m44-m13*m34*m42-m14*m32*m43; 
		float cof13= m12*m23*m44+m13*m24*m42+m14*m22*m43-m12*m24*m43-m13*m22*m44-m14*m23*m42; 
		float cof14= m12*m24*m33+m13*m22*m34+m14*m23*m32-m12*m23*m34-m13*m24*m32-m14*m22*m33;
		
		float cof21= m21*m34*m43+m23*m31*m44+m24*m33*m41-m21*m33*m44-m23*m34*m41-m24*m31*m43;
		float cof22= m11*m33*m44+m13*m34*m41+m14*m31*m43-m11*m34*m43-m13*m31*m44-m14*m33*m41;
		float cof23= m11*m24*m43+m13*m21*m44+m14*m23*m41-m11*m23*m44-m13*m24*m41-m14*m21*m43;
		float cof24= m11*m23*m34+m13*m24*m31+m14*m21*m33-m11*m24*m33-m13*m21*m34-m14*m23*m31;
		
		float cof31= m21*m32*m44+m22*m34*m41+m24*m31*m42-m21*m34*m42-m22*m31*m44-m24*m32*m41;
		float cof32= m11*m34*m42+m12*m31*m44+m14*m32*m41-m11*m32*m44-m12*m34*m41-m14*m31*m42;
		float cof33= m11*m22*m44+m12*m24*m41+m14*m21*m42-m11*m24*m42-m12*m21*m44-m14*m22*m41;
		float cof34= m11*m24*m32+m12*m21*m34+m14*m22*m31-m11*m22*m34-m12*m24*m31-m14*m21*m32;
		
		float cof41= m21*m33*m42+m22*m31*m43+m23*m32*m41-m21*m32*m43-m22*m33*m41-m23*m31*m42;
		float cof42= m11*m32*m43+m12*m33*m41+m13*m31*m42-m11*m33*m42-m12*m31*m43-m13*m32*m41;
		float cof43= m11*m23*m42+m12*m21*m43+m13*m22*m41-m11*m22*m43-m12*m23*m41-m13*m21*m42;
		float cof44= m11*m22*m33+m12*m23*m31+m13*m21*m32-m11*m23*m32-m12*m21*m33-m13*m22*m31;

		return mat4x4(cof11*invDet, cof12*invDet, cof13*invDet, cof14*invDet,
			          cof21*invDet, cof22*invDet, cof23*invDet, cof24*invDet,
					  cof31*invDet, cof32*invDet, cof33*invDet, cof34*invDet,
					  cof41*invDet, cof42*invDet, cof43*invDet, cof44*invDet);
		
	}
	
	float determinant( ) const {
		
		// calculate sub determinants
		return  (m34*m43*m12*m21) - (m33*m44*m12*m21) -
				(m34*m42*m13*m21) + (m32*m44*m13*m21) +
				(m33*m42*m14*m21) - (m32*m43*m14*m21) -
				(m11*m34*m43*m22) + (m11*m33*m44*m22) +
				(m34*m41*m13*m22) - (m33*m41*m14*m22) +
				(m11*m34*m42*m23) - (m11*m32*m44*m23) -
				(m34*m41*m12*m23) + (m32*m41*m14*m23) -
				(m11*m33*m42*m24) + (m11*m32*m43*m24) +
				(m33*m41*m12*m24) - (m32*m41*m13*m24) -
				(m44*m13*m22*m31) + (m43*m14*m22*m31) +
				(m44*m12*m23*m31) - (m42*m14*m23*m31) -
				(m43*m12*m24*m31) + (m42*m13*m24*m31);
    	
	}
	
	// Transformations
	void makeTranslateMatrix(float tx, float ty, float tz) {
		identity();
		m14 = tx; m24 = ty; m34 = tz;
	}
	
	void makeTranslateMatrix(const vec3d &t) {
		makeTranslateMatrix(t.x, t.y, t.z);
	}

	void makeRotXMatrix(float angle) {
		float c = cos(ML_DegToRad(angle)), s = sin(ML_DegToRad(angle));
		identity();
		m22 = c; m23 = -s;
		m32 = s; m33 = c;
	}
	
	void makeRotYMatrix(float angle) {
		float c = cos(ML_DegToRad(angle)), s = sin(ML_DegToRad(angle));
		identity();
		m11 = c;  m13 = s;
		m31 = -s; m33 = c;
	}
	
	void makeRotZMatrix(float angle) {
		float c = cos(ML_DegToRad(angle)), s = sin(ML_DegToRad(angle));
		identity();
		m11 = c; m12 = -s;
		m21 = s; m22 = c;
	}
	
	void makeRotAxis(float angle, const vec3d &v) {
		
		float s=sin(ML_DegToRad(angle)), c=cos(ML_DegToRad(angle)), onemc=1.0f-c;
		float xy, xz, yz;
		vec3d axis(v);

		// guarantee unit rot. axis
		if(1.0f != axis.length()) axis.normalize();

		xy = axis.x*axis.y;
		xz = axis.x*axis.z;
		yz = axis.y*axis.z;
		m11 = axis.x*axis.x*onemc + c;
		m12 = xy*onemc - axis.z*s;
		m13 = xz*onemc + axis.y*s;
		m21 = xy*onemc + axis.z*s;
		m22 = axis.y*axis.y*onemc + c;
		m23 = yz*onemc - axis.x*s;
		m31 = xz*onemc - axis.y*s;
		m32 = yz*onemc + axis.x*s;
		m33 = axis.z*axis.z*onemc + c;
		m41= m42= m43= m14= m24= m34= 0.0f;
		m44= 1.0f;
	}
	
	void makeScaleMatrix( float sx, float sy, float sz) {
		identity();
		m11 = sx; m22 = sy; m33 = sz;
	}

	void makePersProjMatrixFov( float fovy, float aspect,
							    float near, float far) {
		
		float tanF = tan(0.5f*ML_DegToRad(fovy));
        assert(tanF>ML_EPSILON);
		float cotF = 1.0f/tanF;
		float deltaZ= near-far;

		identity();
		m11 = cotF/aspect;
		m22 = cotF;
		m33 = (far+near)/deltaZ;
		m34=  (2.0f*far*near)/deltaZ;
		m43 = -1.0f;
		m44 = 0.0f;

	}

	void makePersProjMatrix( float left,   float right,
		                     float bottom, float top,
							 float near,   float far) {

		float deltaX, deltaY, deltaZ;
		deltaX=right-left;
		deltaY=top-bottom;
		deltaZ=far-near;

		identity();
		m11= 2.0f*near/deltaX;
		m13= (right+left)/deltaX;
		m22= 2.0f*near/deltaY;
		m23= (top+bottom)/deltaY;
		m33= -(far+near)/deltaZ;
		m34= -(2.0f*far*near)/deltaZ;
		m43= -1.0f;
	
	}

	void makeOrthoProjMatrix( float left,   float right,
		                      float bottom, float top,
							  float near,   float far ) {
		float deltaX, deltaY, deltaZ;
		deltaX=right-left;
		deltaY=top-bottom;
		deltaZ=far-near;

		identity();
		m11= 2.0f/deltaX;
		m14= -(right+left)/deltaX;
		m22= 2.0f/deltaY;
		m24= -(top+bottom)/deltaY;
		m33= -2.0f/deltaZ;
		m34 = -(far+near)/deltaZ;

	}
	
	void makeOrtho2DProjMatrix( float left,  float right,
								float bottom, float top ) {
		
		makeOrthoProjMatrix(left, right, bottom, top, -1.0f, 1.0f);
	}
	
	void makeLookAtRH( float eyeX,  float eyeY, float eyeZ,
		               float lookX, float lookY, float lookZ,
					   float upX,   float upY,  float upZ) {
		
		vec3d eye(eyeX, eyeY, eyeZ);
		vec3d look(lookX, lookY, lookZ);
		vec3d up(upX, upY, upZ);
		vec3d f, u, s;
		
		// f = look - eye
		f = look-eye;
		f.normalize();
		
		if(1.0f != up.length() ) up.normalize();
		
		// s = f x up
		s = f.cross(up);
		
		// u = s x f
		u = s.cross(f);

		m11=  s.x;  m12= s.y;  m13= s.z;  m14= -eye.dot(s);
		m21=  u.x;  m12= u.y;  m13= u.z;  m14= -eye.dot(u);
		m11= -f.x;  m12= -f.y; m13= -f.z; m14= -eye.dot(f);
		m11=  0.0f; m12= 0.0f; m13= 0.0f; m14= 1.0f;
	}

	void makeLookAtRH(const vec3d &eye, const vec3d &look, const vec3d &up) {
		makeLookAtRH( eye.x, eye.y, eye.z, 
                      look.x, look.y, look.z,
					  up.x, up.y, up.z );
	}

	// Query methods
	bool isSymmetric( ) const { return *this == getTranspose(); }
	
	bool isDiagonal ( ) const {
		return (  m12==0.0f && m13==0.0f && m14==0.0f &&
			      m21==0.0f && m23==0.0f && m24==0.0f &&
			      m31==0.0f && m32==0.0f && m34==0.0f &&
			      m41==0.0f && m42==0.0f && m43==0.0f );
	}
	
	bool isIdentity ( ) const {
		return ( m11==1.0f && m22==1.0f && m33==1.0f && m44==1.0f) &&
	           ( m12==0.0f && m13==0.0f && m14==0.0f) &&
			   ( m21==0.0f && m23==0.0f && m24==0.0f) &&
			   ( m31==0.0f && m32==0.0f && m34==0.0f) &&
			   ( m41==0.0f && m42==0.0f && m43==0.0f); 
	}

	union {
		struct {
			float m11, m12, m13, m14;
			float m21, m22, m23, m24;
			float m31, m32, m33, m34;
			float m41, m42, m43, m44;
		};
		float m[16];
	};
};


//===========================================================================
//
// Quaternion Class
//
//===========================================================================

//struct quat {
//	// Constructors
//	explicit Quaternion();
//			 Quaternion(float X, float Y, float Z, float W);
//			 Quaternion(const Vector3D &Imaginer, float Real);
//	explicit Quaternion(const Vector4D &QuadVec);
//	explicit Quaternion(const float *quat_array);
//			 Quaternion(const Quaternion &quaternion);
//	
//	// Unary operators
//	Quaternion  operator + ( ) const;
//	Quaternion  operator - ( ) const;
//	
//	// Binary operators
//		   Quaternion operator + (const Quaternion &quaternion) const;
//		   Quaternion operator - (const Quaternion &quaternion) const;
//		   Quaternion operator * (const Quaternion &quaternion) const;
//		   Quaternion operator * (float scalar) const; // Multiply from a scalar to the rhs
//	friend Quaternion operator * (float scalar, const Quaternion &quaternion); // Multiply from a scalar to the lhs
//		   Quaternion operator / (float scalar) const;
//		
//	// Assignment operators
//	Quaternion& operator += (const Quaternion& quaternion);
//	Quaternion& operator -= (const Quaternion& quaternion);
//	Quaternion& operator *= (const Quaternion &quaternion);
//	Quaternion& operator *= (float scalar);
//	Quaternion& operator /= (float scalar);
//	
//	// Casting
//	operator float* ();
//	operator const float* ();
//	
//	// Equality operators
//	bool operator == (const Quaternion &quaternion) const;
//	bool operator != (const Quaternion &quaternion) const;
//	void set(float X, float Y, float Z, float W);
//	void SetIdentityQuaternion();
//	
//	// Quaternion consists of 3 Imaginer and 1 Real number.
//	// Q = x*i + y*j + z*k + real or Q = [img, real] where img = <x, y, z>  
//	union {
//		struct {
//			float x;
//			float y;
//			float z;
//			float real;
//		};
//		float data[4];
//	};
//
//};


//===========================================================================
//
// Ray Class
//
//===========================================================================
struct ray {
	// Constructors
	explicit ray() {/* Empty constructor is preferred due to perf. reasons */ }
			 ray( float ox, float oy, float oz,
				  float dx, float dy, float dz ) : vOrg(ox, oy, oz),
				                                   vDir(dx, dy, dz) { }
			
			ray(vec3d &org, vec3d &dir) : vOrg(org), vDir(dir) { }
			ray(const ray &r) : vOrg(r.vOrg), vDir(r.vDir) { }
	
	// Accessors
	void set(const vec3d &org, const vec3d &dir) {
		 vOrg = org; vDir=dir;
	}
	
	void setOrigin(vec3d &org) {
		vOrg = org;
	}
	
	void  setDirection(vec3d &dir) {
		vDir = dir;
		// always normalized
		dir.normalize();
	}

	vec3d getOrigin( ) const { return vOrg; }
	vec3d getDir   ( ) const { return vDir; }

	// Ray is represented by the following equation r(t) = vOrg + t*vDir
	// Here vOrg is the ray postion and vDir is the ray's direction vector.
	vec3d  vOrg;
	vec3d  vDir;
	
};


//===========================================================================
//
// Plane Class
//
//===========================================================================
struct plane {
	// Constructors
	explicit plane() {/* Empty constructor is preferred due to perf. reasons */ }
			 plane(float nx, float ny, float nz, float dist) {
				 normal.set(nx, ny, nz);
				 d= dist;
	         }
			 // construct plane from the given three points
			 plane(const vec3d &p1, const vec3d &p2, const vec3d &p3) {
				 normal.x= p1.y*(p2.z-p3.z)+p2.y*(p3.z-p1.z)+p3.y*(p1.z-p2.z);
				 normal.y= p1.z*(p2.x-p3.x)+p2.z*(p3.x-p1.x)+p3.z*(p1.x-p2.x);
				 normal.z= p1.x*(p2.y-p3.y)+p2.x*(p3.y-p1.y)+p3.x*(p1.y-p2.y);
				 d= -(p1.x*(p2.y*p3.z-p3.y*p2.z)+p2.x*(p3.y*p1.z-p1.y*p3.z)+ p3.x*(p1.y*p2.z-p2.y*p1.z));
			 }
	explicit plane(const float *parray) {
				assert(parray);
				normal.set(parray[0], parray[1], parray[2]);
				d= parray[3]; 
			 }
			 plane(const plane &pln) {
				 normal.set(pln.normal.x, pln.normal.y, pln.normal.z);
				 d= pln.d;
			 }
			 
			 ~plane( ) { }
	
	// Accessors
	vec3d getNormal()const { return normal; }
	float getDistance()const { return d; }
	vec4d getCoeffs()const { return vec4d(normal, d); }
	void  setNormal(vec3d &n) { normal=n; }
	void  setDistance(float dist) { d = dist; }
	int   pointPlaneRelation(const vec3d &p) {
		// calculate dot(normal, p) + d and determine its sign
		float s = p.dot(normal) + d; 
		return s > 0.0f ? PPR_FRONT : s==0.0f ? PPR_ON : PPR_BACK;  
	}
	
	// Point-Plane relationship type
	enum {
		PPR_FRONT=0,
		PPR_BACK,
		PPR_ON
	};

	// A Plane is represented by the equation: Ax+By+Cz+d=0
	// where normal=[A B C] .
	vec3d  normal;
	float  d;

};


//===========================================================================
//
// AABB class
//
//===========================================================================
struct AABB {
	explicit AABB() { min.setToZeroVector(); max.setToZeroVector(); }
	         AABB(const vec3d &minp, const vec3d &maxp) {
		          min=minp; max=maxp;
	         }

	vec3d  getCenter()       { return (min + max)*0.5f; }
	vec3d  getSize()         { return  max - min;       }
	float  getWidth()  const { return max.x - min.x;    }
	float  getHeight() const { return max.y - min.y;    }
	float  getLength() const { return max.z - min.z;    }
	void   setEmpty() {
		   min.setToZeroVector();
		   max.setToZeroVector();
	}

	enum PLANE_RELATION {
		PR_INSIDE=0,
		PR_OUTSIDE,
		PR_INTERSECT
	};
		
	// plane intersection relation
	// adapted from Real Time Rendering,3rd.ed. p.756
	PLANE_RELATION planeIntersect(const plane &p) {
		vec3d c = (max + min)*0.5f;
		vec3d h = (max - min)*0.5f;
		float e = h.x*fabs(p.normal.x) +
			      h.y*fabs(p.normal.y) +
				  h.z*fabs(p.normal.z);

		float s = c.dot(p.normal) + p.d;

		if(s - e > 0) return PR_OUTSIDE;
		if(s + e < 0) return PR_INSIDE;

		return PR_INTERSECT;
	}

    
	// extreme points of the box
	vec3d  min;
	vec3d  max;
};


//===========================================================================
// 2D Bounding Rectangle
//
// (x,y)............(x+w,y)
//   .					.
//   .					.
//	 .					.
// (x,y+h).........(x+w,y+h)
// 
//===========================================================================
struct rectangle {
	   explicit rectangle( ) { }
	            rectangle(float x, float y, float w, float h) :
	                      pos(x,y) , width(w), height(h) { }
			    rectangle(const vec2d &position, const vec2d &dimension) :
						  pos(position), width(dimension.x), height(dimension.y) { }
	            rectangle(const rectangle &cpy) {
					pos=cpy.pos;
					width=cpy.width;
					height=cpy.height;
			    }

				rectangle& operator = (const rectangle &cpy) {
					pos=cpy.pos; width=cpy.width; height=cpy.height;
					return *this;
				}

	// Accessors
   float  getLeft( )   const { return pos.x;                            }
   float  getRight( )  const { return pos.x+width;                      }
   float  getTop( )    const { return pos.y;                            }
   float  getBottom( ) const { return pos.y+height;                     }
   vec2d  getCenter( ) const { return vec2d(width/2.0f, height/2.0f);   }
   void   setEmpty( )        { pos.setToZeroVector(); width=height = 0; }

    // top left corner position
	vec2d   pos;   
    // dimensions
	float   width; 
	float   height;
};

//===========================================================================
//
// Frustum 
//
//===========================================================================
struct frustum {
	   explicit frustum () { }
	   frustum(const mat4x4 matModelView, const mat4x4 matProj, bool normalize) {
	   
		  extractPlanes(matModelView*matProj, normalize);
	   }
	   
	   frustum(const mat4x4 &matModelViewProj, bool normalize) {
		   extractPlanes(matModelViewProj, normalize);
	   }

	   frustum(const frustum &copy) {
			for(int i=0;i<6;++i) 
				planes[i] = copy.planes[i];
	   }


	   // Adapted from the paper: "Fast Extraction of Viewing
	   // Frustum Planes from the World-View-Projection Matrix"
	   // Note: It assumes the matrix 'mat' is a projection matrix.
	   //       In the case of a concataneted matrix(ModelView*Proj) as 
	   //       in OpenGL, the extracted planes will be in the 
	   //       object space instead of eye space.
	   void extractPlanes(const mat4x4 &mat, bool normalizeplanes) {
			
		   // left plane
		   plane[0] = plane(mat.m41+mat.m11,
			                mat.m42+mat.m12,
							mat.m43+mat.m13,
							mat.m44+mat.m14 );

		   // right plane
		   plane[1] = plane(mat.m41-mat.m11,
			                mat.m42-mat.m12,
							mat.m43-mat.m13,
							mat.m44-mat.m14 );
		
		   // bottom plane
		   plane[2] = plane(mat.m41+mat.m21,
			                mat.m42+mat.m22,
							mat.m43+mat.m23,
							mat.m44+mat.m24 );
		
		   // top plane
		   plane[3] = plane(mat.m41-mat.m21,
			                mat.m42-mat.m22,
							mat.m43-mat.m23,
							mat.m44-mat.m24 );
		
		   // near plane
		   plane[4] = plane(mat.m41+mat.m31,
			                mat.m42+mat.m32,
							mat.m43+mat.m33,
							mat.m44+mat.m34 );
		
		   // far plane
		   plane[5] = plane(mat.m41-mat.m31,
			                mat.m42-mat.m32,
							mat.m43-mat.m33,
							mat.m44-mat.m34 );
		
		   // normalize if necessary
		   if(normalizeplanes) {
			   for(int k=0;k<6;k++)
				   plane[k].normal.normalize();
		   }

	   }

	   void extractPlanes(const mat4x4 &viewmodel, const mat4x4 &proj,
		                  bool normalize) {
		   extractPlanes(proj*viewmodel, normalize);
	   }
	   
	   
	   // frustum-AABB relation status
	   enum AABB_Relation {
		   AABB_INSIDE=0,
		   AABB_OUTSIDE,
		   AABB_INTERSECT
	   };


	   // Frustum-AABB intersection method
	   // Adapted from RTR,3ed, p.777
 	   AABB_Relation  aabbIntersect(const AABB &aabb) {
		   AABB::PLANE_RELATION  result;
		   bool intersecting = false;
		   
			// test all 6 planes with AABB	
		   for(int i=0;i<6;++i) {
			   result = aabb.planeIntersect(plane[i]);
	   		   if(result == PR_OUTSIDE) return AABB_OUTSIDE;
				else if(result == PR_INTERSECT)
						intersecting=true;
			   if(intersecting) return AABB_INTERSECT;
			   else return AABB_INSIDE;
		   }

	   }
  
		// six plane names 
	   enum {
		   LEFT_PLANE=0,
		   RIGHT_PLANE,
		   BOTTOM_PLANE,
		   TOP_PLANE,
		   NEAR_PLANE,
		   FAR_PLANE,
	   };

		
	   // six frustum planes
	   plane planes[6];
};

//===========================================================================
//
// RGB Floating Point Color Class
//
//===========================================================================
struct colorRGB {
	explicit colorRGB() {/* Empty constructor is preferred due to perf. reasons */ }
			 colorRGB(float Red, float Green, float Blue) : R(Red), G(Green), B(Blue) { }
			 explicit colorRGB(const float *parray) {assert(parray); R=parray[0]; G=parray[1]; B=parray[2]; }
			 explicit colorRGB(const vec3d &colorvec) { R= colorvec.x; G= colorvec.y; B= colorvec.z; }
			 colorRGB(const colorRGB &copy) : R(copy.R), G(copy.G), B(copy.B) { }
	//explicit ColorRGB(const unsigned long  value);
		
	// binary operators
	colorRGB operator + (const colorRGB &c) const {
		return colorRGB(R+c.R, G+c.G, B+c.B);
	}

	colorRGB operator - (const colorRGB &c) const {
		return colorRGB(R+c.R, G+c.G, B+c.B);
	}
    
	colorRGB operator * (float scalar) const {
		return colorRGB(R*scalar, G*scalar, B*scalar);
	}

	colorRGB operator * (const colorRGB &c) const {
		return colorRGB(R*c.R, G*c.G, B*c.B);
	}

	colorRGB operator / ( float scalar) const {
		assert(scalar>ML_EPSILON);
		return colorRGB(R/scalar, G/scalar, B/scalar);
	}

	// unary operators
	colorRGB operator + () {
		return colorRGB(R++, G++, B++);
	}

	colorRGB operator - () {
		return colorRGB(R--, G--, B--);
	}

	// Casting operators
	operator float* () { return (float*)&R; }
	operator const float* ( ) const { return (const float*)&R; }        
		    
    // Assignment operators
	colorRGB& operator += (const colorRGB &c ) {
		return *this = *this + c;
	}
    
	colorRGB& operator -= ( const colorRGB  &c ) {
		return *this = *this - c;
	}
    
	colorRGB& operator *= ( float scalar ) {
		return *this = *this * scalar;
		
	}

	colorRGB& operator *= ( const colorRGB &c ) {
		return *this = *this * c;
	}

    
	colorRGB& operator /= ( float scalar ) {
		assert(scalar>ML_EPSILON);
		return *this = *this / scalar;
	}

	colorRGB& operator = (const colorRGB &c) {
		R= c.R; G= c.G; B= c.B;
		return *this;
	}

	friend colorRGB operator * ( float scalar, const colorRGB &c ) {
		return colorRGB(scalar*c.R, scalar*c.G, scalar*c.B);
	}

	void set(float _R, float _G, float _B) {
		R= _R; G= _G; B= _B;
	}
	
	void set(const colorRGB &color) {
		R= color.R; G= color.G; B= color.B;
	}

	// Equility 
	bool operator == (const colorRGB &c) const {return ( R == c.R && G == c.G && B == c.B ); }
	bool operator != (const colorRGB &c) const {return ( R != c.R || G != c.G || B != c.B ); }
	bool operator >= (const colorRGB &c) const {return ( R >= c.R && G >= c.G && B >= c.B ); }
	bool operator <= (const colorRGB &c) const {return ( R <= c.R && G <= c.G && B <= c.B ); }
	bool operator >  (const colorRGB &c) const {return ( R  > c.R && G  > c.G && B >  c.B ); }
	bool operator <  (const colorRGB &c) const {return ( R  < c.R && G  < c.G && B <  c.B ); } 
	
	// FP Color components : Red, Green Blue in the range [0.0, 1.0].
	union {
		struct {
			float R;
			float G;
			float B;
		};
		float color[3];
	};
};





//===========================================================================
//
// Standart Basis Vectors
//
//===========================================================================
static const vec2d    ZeroVec2D(0.0f, 0.0f);
static const vec2d    UNIT_XVec2D(1.0f, 0.0f);
static const vec2d    UNIT_YVec2D(0.0f, 1.0f);
static const vec2d    NEGATIVE_UNIT_XVec2D(-1.0f, 0.0f);
static const vec2d    NEGATIVE_UNIT_YVec2D(0.0f, -1.0f);
static const vec2d    UNIT_SCALEVec2D(1.0f, 1.0f);
static const vec3d    ZeroVec3D(0.0f, 0.0f, 0.0f);
static const vec3d    UNIT_XVec3D(1.0f, 0.0f, 0.0f);
static const vec3d    UNIT_YVec3D(0.0f, 1.0f, 0.0f);
static const vec3d    UNIT_ZVec3D(0.0f, 0.0f, 1.0f);
static const vec3d    NEGATIVE_UNIT_XVec3D(-1.0f, 0.0f, 0.0f);
static const vec3d    NEGATIVE_UNIT_YVec3D(0.0f, -1.0f, 0.0f);
static const vec3d    NEGATIVE_UNIT_ZVec3D(0.0f, 0.0f, -1.0f);
static const vec3d    UNIT_SCALEVec3D(1.0f, 1.0f, 1.0f);
static const vec4d    ZeroVec4D(0.0f, 0.0f, 0.0f, 0.0f);
static const vec4d    UNIT_XVec4D(1.0f, 0.0f, 0.0f, 0.0f);
static const vec4d    UNIT_YVec4D(0.0f, 1.0f, 0.0f, 0.0f);
static const vec4d    UNIT_ZVec4D(0.0f, 0.0f, 1.0f, 0.0f);
static const vec4d    NEGATIVE_UNIT_XVec4D(-1.0f, 0.0f, 0.0f, 0.0f);
static const vec4d    NEGATIVE_UNIT_YVec4D(0.0f, -1.0f, 0.0f, 0.0f);
static const vec4d	  NEGATIVE_UNIT_ZVec4D(0.0f, 0.0f, -1.0f, 0.0f);
static const vec4d	  UNIT_SCALEVec4D(1.0f, 1.0f, 1.0f, 0.0f);
static const plane    POSX_PLANE(0.0f,  1.0f, 0.0f, 0.0f); 
static const plane    NEGX_PLANE(0.0f, -1.0f, 0.0f, 0.0f); 
static const plane    POSY_PLANE(1.0f,  0.0f, 0.0f, 0.0f); 
static const plane    NEGY_PLANE(0.0f, -1.0f, 0.0f, 0.0f); 


//===========================================================================
//
// Primary colors
//
//===========================================================================
static const colorRGB BLACK_COLOR(0.0f, 0.0f, 0.0f);
static const colorRGB WHITE_COLOR(1.0f, 1.0f, 1.0f);
static const colorRGB RED_COLOR(1.0f, 0.0f, 0.0f);
static const colorRGB BLUE_BLACK_COLOR(0.0f, 0.0f, 1.0f);
static const colorRGB GREEN_COLOR(0.0f, 1.0f, 0.0f);
static const colorRGB YELLOW_COLOR(1.0f, 1.0f, 0.0f);
static const colorRGB MAGENTA_COLOR(1.0f, 0.0f, 1.0f);

} // end of namespace MATHLINEAR

#endif // end of "MathLinear.h"
