#ifndef LOGREG_HPP
#define LOGREG_HPP

#include <math.h>
#include <vector>
#include <stdio.h>

class Logreg
{
	public:
		Logreg();
		Logreg(double period_);
		Logreg(double offset_, double scale_);

		double train();
		void add_measurement(double time, bool positive, float weight);
		double evaluate(double time) const;

		float my_weight;

		void save(FILE* file, bool lossy = false);
		void load(FILE* file);
		void exportToArray(double* array, int maxLen, int& pos);
		void importFromArray(double* array, int len, int& pos);

		void print();

	private:

		class Vec3D {
			public:
				double x,y,z;
				Vec3D(double x_, double y_, double z_) :
					x(x_),
					y(y_),
					z(z_)
				{

				}

				Vec3D operator *(double s) const {
					return Vec3D(x*s, y*s, z*s);
				}

				double operator *(Vec3D v) const {
					return x*v.x + y*v.y + z*v.z;
				}

				Vec3D operator +(Vec3D v) const {
					return Vec3D(x+v.x, y+v.y, z+v.z);
				}
				Vec3D operator +=(Vec3D v) {
					x += v.x;
					y += v.y;
					z += v.z;
					return *this;
				}
				Vec3D operator -(Vec3D v) const {
					return Vec3D(x-v.x, y-v.y, z-v.z);
				}
				Vec3D operator -=(Vec3D v) {
					x -= v.x;
					y -= v.y;
					z -= v.z;
					return *this;
				}

				double norm() const {
					return sqrt(x*x + y*y + z*z);
				}
		};

		Vec3D w;

		class Measurement {
			public:
				Measurement(Vec3D x_, double y_, double weight_) :
					x(x_),
					y(y_),
					weight(weight_)
				{

				}
				Vec3D x;
				double y;
				double weight;
		};

		std::vector<Measurement> points;
		double total_weight;
		double period;
		double scale;
		double offset;
		bool parabolic;
		bool unlimited;

		double loss(Vec3D at) const;
		Vec3D gradient(Vec3D at) const;
		Vec3D time_to_point(double time) const;
};

#endif // LOGREG_HPP
