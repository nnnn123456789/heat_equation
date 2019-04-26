
#include <array>
#include <iostream>
#include <cmath>
#include <memory>

typedef double _Data;

constexpr long long xN = 1000;
constexpr double xstep = 1.0 / xN;
constexpr long long tN = 1000;
constexpr long long tcount = 10000;
constexpr double tstep = 1.0 / xN;
constexpr double pi = 3.141592653589793;

template<int n>
std::array<_Data, n> solve_catch(const std::array<_Data, n - 1> & _d, const std::array<_Data, n>& _c, const std::array<_Data, n - 1> & _u, const std::array<_Data, n >& _b)
{
	auto pd = new std::array<_Data, n - 1>;
	auto pc = new std::array<_Data, n>;
	auto pu = new std::array<_Data, n - 1>;
	auto pb = new std::array<_Data, n>;
	auto& d = *pd; d = _d;
	auto& c = *pc; c = _c;
	auto& u = *pu; u = _u;
	auto& b = *pb; b = _b;
	std::array<_Data, n> x;
	for (long long i = 1; i < n; i++)
	{
		d[i - 1] /= c[i - 1];
		c[i] -= u[i - 1] * d[i - 1];
		b[i] -= b[i - 1] * d[i - 1];
	}
	x[n - 1] = b[n - 1] / c[n - 1];
for (long long i = n - 2; i >= 0; i--)
{
	x[i] = (b[i] - u[i] * x[i + 1]) / c[i];
}
delete pb, pc, pu, pd;
return x;
}


int main1([[maybe_unused]] int argc = 1, [[maybe_unused]] char* argv[] = nullptr)
{
	constexpr double h = xstep;
	constexpr int J = 10;
	constexpr int k = 50;
	std::unique_ptr<std::array<double, xN + 1 >> plast(new std::array<double, xN + 1>);
	std::unique_ptr<std::array<double, xN + 1>> pcurrent(new std::array<double, xN + 1>);
	std::array<double, xN + 1> & last = *plast;
	std::array<double, xN + 1> & current = *pcurrent;
	std::array<double, xN / J - 1> c;
	std::array<double, xN / J - 2> u, d;
	constexpr double M = -h * h / tstep;
	for (int i = 0; i < (xN / J - 2); i++)
	{
		u[i] = d[i] = 1;
		c[i] = M - 2;
	}
	c[(xN / J - 2)] = M - 2;
	for (int i = 0; i <= xN; i++)
	{
		last[i] = std::sin(pi * i * h);
	}
	for (int i = 0; i < tcount; i++)
	{
		for (long long j = xN / J; j <= xN - xN / J; j += xN / J)
		{
			double RHS = (last[j - k] - 2 * last[j] + last[j + k]) / (k * h * k * h);
			current[j] = last[j] + tstep * RHS;
		}
		current[0] = current[xN] = 0;
		for (long long j = 0; j < J; j++)
		{
			std::array<double, xN / J - 1> right;
			for (long long l = 1; l <= xN / J - 1; l++)
				right[l - 1] = M * last[xN / J * j + l];
			right[0] -= current[j * xN / J];
			right[xN / J - 2] -= current[(j + 1) * xN / J];
			auto cur = solve_catch<xN / J - 1>(d, c, u, right);
			for (long long l = 1; l <= xN / J - 1; l++)
				current[xN / J * j + l] = cur[l - 1];
		}
		if (0 == i % 200)
		{
			printf("t=%f\n", i * tstep);
			for (int k = 0; k <= xN; k += 50)
			{
				printf("%.10f\t", last[k]);
				if (k % 250 == 0)
					printf("\n");
			}
			printf("\n\n");
		}
		last = current;

	}
	system("pause");
	return 0;
}

int main2([[maybe_unused]] int argc = 1, [[maybe_unused]] char* argv[] = nullptr)
{
	constexpr double h = xstep;
	constexpr int k = 50;
	std::unique_ptr<std::array<double, xN + 1>> plast(new std::array<double, xN + 1>);
	std::unique_ptr<std::array<double, xN + 1>> pcurrent(new std::array<double, xN + 1>);
	std::array<double, xN + 1> & last = *plast;
	std::array<double, xN + 1> & current = *pcurrent;
	std::unique_ptr<std::array<double, xN - 1>> pc(new std::array<double, xN - 1>);
	std::unique_ptr<std::array<double, xN - 2>> pu(new std::array<double, xN - 2>);
	std::unique_ptr<std::array<double, xN - 2>> pd(new std::array<double, xN - 2>);
	auto & c = *pc;
	auto & u = *pu;
	auto & d = *pd;
	constexpr double M = -h * h / tstep;
	for (int i = 0; i < xN - 2; i++)
	{
		u[i] = d[i] = 1;
		c[i] = M - 2;
	}
	c[xN - 2] = M - 2;
	for (int i = 0; i <= xN; i++)
	{
		last[i] = std::sin(pi * i * h);
	}
	for (int i = 0; i < tcount; i++)
	{

		current[0] = current[xN] = 0;
		std::unique_ptr<std::array<double, xN - 1>> pright(new std::array<double, xN - 1>);
		auto & right = *pright;
		for (long long l = 1; l <= xN - 1; l++)
			right[l - 1LL] = M * last[l];
		auto cur = solve_catch<xN - 1>(d, c, u, right);
		for (long long l = 1; l <= xN - 1; l++)
			current[l] = cur[l - 1LL];

		if (0 == i % 200)
		{
			printf("t=%f\n", i * tstep);
			for (int k = 0; k <= xN; k += 50)
			{
				printf("%.10f\t", last[k]);
				if (k % 250 == 0)
					printf("\n");
			}
			printf("\n\n");
		}
		last = current;
	}
	system("pause");
	return 0;
}

int main3([[maybe_unused]] int argc = 1, [[maybe_unused]] char* argv[] = nullptr)
{
	auto u = [](double t, double x) noexcept ->double
	{
		return exp(-pi * pi * t) * sin(pi * x);
	};
	for (int i = 0; i <= tcount; i++)
	{
		if (0 == i % 200)
		{
			printf("t=%f\n", i * tstep);
			for (int k = 0; k <= xN; k += 50)
			{
				printf("%.10f\t", u(i * tstep, k * xstep));
				if (k % 250 == 0)
					printf("\n");
			}
			printf("\n\n");
		}
	}
	system("pause");
	return 0;
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[])
{
	main1(argc, argv);
	main2(argc, argv);
	main3(argc, argv);
	return 0;
}
