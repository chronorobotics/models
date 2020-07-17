#ifndef SQUAREERROR_H
#define SQUAREERROR_H

#include <vector>

class SquareError
{
	public:
		SquareError(int size_);

		std::vector<double> forward(const std::vector<double>& x, const std::vector<double>& y, int batch_size);
		std::vector<double> backward();

		void set_size(int new_size);

	private:
		int size;
		std::vector<double> last_x;
		std::vector<double> last_y;
		int bs;
};

#endif // SQUAREERROR_H
