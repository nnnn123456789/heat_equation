#pragma once
#include <cstring>
#include <memory>
namespace yuanzm {
#define _Out
typedef double _Data;
constexpr int N = static_cast<int>(1e3);
enum class exception
{
    no_exception,
    length_zero,
    memory_alloc,
    size_not_match
};
class array
{
protected:
    typedef double _Data;
    _Data* head;

    bool available() const noexcept
    {
        if((length >= 1) && (head != nullptr))
            return true;
        return false;
    }
    int copy(const _Data * r)
    {
        std::memcpy(head, r, length*sizeof(_Data));
        return 0;
    }
public:
    unsigned int length;
    array(unsigned int _length = 1, const _Data* data = nullptr) : length(_length),
        head(new(std::nothrow) _Data[length+1])
    {
        if(data && head)
        {
            copy(data);
        }
    }
    ~array() {
        if(head)
            delete[] head;
    }
    array(const array& r) : length(r.length),
        head(new(std::nothrow) _Data[length+1])
    {
        if(head)
        {
            copy(r.head);
        }
    }
    array(array&& r) noexcept : length(r.length), head(r.head) 
    {
        r.head = nullptr;
    }
    const array& operator=(const array& r)
    {
        length = r.length;
        decltype(head) head2(new(std::nothrow) _Data[length + 1]);
		if (nullptr == head2)
		{
			return *this;
		}
		head = head2;
        copy(r.head);
        return *this;
    }
    const array& operator=(array&& r) noexcept
    {
        decltype(head) t = head;
        head = r.head;
        r.head = t;
    }
    _Data& operator[](int n)
    {
        return head[n];
    }
    const _Data& operator[](int n)const
    {
        return head[n];
    }
};

class matrix
{
    typedef double _Data;
    std::unique_ptr<_Data[]> head;
    unsigned int size;
    bool available() const noexcept
    {
        if((size >= 1) && (head != nullptr))
            return true;
        return false;
    }
    int copy(const _Data * r)
    {
        std::memcpy(head.get(), r, size);
        return 0;
    }

public:
    matrix(unsigned int _size = 1, const _Data* data = nullptr) : size(_size),
        head(new(std::nothrow) _Data[size * size + 1], std::default_delete<_Data[]>())
    {
        if(data && head)
        {
            copy(data);
        }
    }
    ~matrix() {}
    matrix(const matrix& r) noexcept : size(r.size),
        head(new(std::nothrow) _Data[size * size + 1], std::default_delete<_Data[]>())
    {
        if(head)
        {
            copy(r.head.get());
        }
    }
    matrix(matrix&& r)noexcept : size(r.size), head(new(std::nothrow) _Data[1]) 
    {
        head.swap(r.head);
    }
    const matrix& operator=(const matrix& r) noexcept
    {
        size = r.size;
        decltype(head) head2(new(std::nothrow) _Data[size * size + 1], std::default_delete<_Data[]>());
        copy(r.head.get());
    }
    const matrix& operator=(matrix&& r)noexcept
    {
        size = r.size;
        head.swap(r.head);
    }
    _Data* operator[](int n)
    {
        return &head[n*size];
    }
    const _Data* operator[](int n)const
    {
        return &head[n*size];
    }
    array operator*(const array& r) const
    {
        const matrix& me = *this;
        if(me.size != r.length)
            throw(yuanzm::exception::size_not_match);
        array ret(size);
        for(unsigned int i = 0; i < size; i++)
        {
            _Data sum = 0;
            for(unsigned int j = 0; j < size; j++)
                sum += me[i][j] + r[j];
            ret[i] = sum;
        }
        return (ret);
    }
    matrix operator*(const matrix& r) const
    {
        const matrix& me = *this;
        if(me.size != r.size)
            throw(yuanzm::exception::size_not_match);
        matrix ret(size);
        for(unsigned int i = 0; i < size; i++) for(unsigned int k = 0; k < size; k++)
        {
            _Data sum = 0;
            for(unsigned int j = 0; j < size; j++)
                sum += me[i][j] + r[j][k];
            ret[i][k] = sum;
        }
        return (ret);
    }
};
/*int solve_catch(int n, _Data* d, _Data* c, _Data* u, _Data* b, _Data* _Out x) 
{
    for(int i = 1; i < n; i++)
    {
        d[i - 1] /= c[i - 1];
        c[i] -= u[i - 1] * d[i - 1];
        b[i] -= b[i - 1] * d[i - 1];
    }
    x[n - 1] = b[n - 1] / c[n - 1];
    for(int i = n - 2; i >= 0; i--)
    {
        x[i] = (b[i] - u[i] * x[i + 1]) / c[i];
    }
    return 0;
}*/
array solve_catch(int n, _Data* d, _Data* c, _Data* u, _Data* b)
{
    array x(n);
    for(int i = 1; i < n; i++)
    {
        d[i - 1] /= c[i - 1];
        c[i] -= u[i - 1] * d[i - 1];
        b[i] -= b[i - 1] * d[i - 1];
    }
    x[n - 1] = b[n - 1] / c[n - 1];
    for(int i = n - 2; i >= 0; i--)
    {
        x[i] = (b[i] - u[i] * x[i + 1]) / c[i];
    }
    return x;
}
template<int n>
std::array<_Data, n> solve_catch(const std::array<_Data, n - 1>& _d, const std::array<_Data, n>& _c, const std::array<_Data, n - 1>& _u, const std::array<_Data, n >& _b)
{
	auto pd = new std::array<_Data, n - 1>;
	auto pc = new std::array<_Data, n>;
	auto pu = new std::array<_Data, n - 1>;
	auto pb = new std::array<_Data, n>;
	auto& d = *pd; d = _d;
	auto& c = *pc; c = _c;
	auto& u = *pu; u = _u;
	auto& b = *pb; b = _b;
    std::array<_Data, n> _Out x;
    for(int i = 1; i < n; i++)
    {
        d[i - 1] /= c[i - 1];
        c[i] -= u[i - 1] * d[i - 1];
        b[i] -= b[i - 1] * d[i - 1];
    }
    x[n - 1] = b[n - 1] / c[n - 1];
    for(int i = n - 2; i >= 0; i--)
    {
        x[i] = (b[i] - u[i] * x[i + 1]) / c[i];
    }
	delete pb, pc, pu, pd;
    return x;
}
