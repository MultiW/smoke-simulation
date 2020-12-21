#include <vector>
#include <stdio.h>

/*
* Based on: https://stackoverflow.com/questions/18062463/what-is-the-most-efficient-way-to-initialize-a-3d-vector
* - Using a contiguous array provides good spatial locality
*/

template <typename T>
class List3d {
public:
    List3d(size_t d1 = 0, size_t d2 = 0, size_t d3 = 0, T const& t = T()) :
        d1(d1), d2(d2), d3(d3), data(d1* d2* d3, t)
    {}

    T& operator()(size_t i, size_t j, size_t k)
    {
        return data[i*d2*d3 + j*d3 + k];
    }

    T const& operator()(size_t i, size_t j, size_t k) const
    {
        return data[i*d2*d3 + j*d3 + k];
    }

    void resize(size_t d1, size_t d2, size_t d3)
    {
        this->d1 = d1;
        this->d2 = d2;
        this->d3 = d3;
        data.resize(d1 * d2 * d3);
    }

    size_t size(int dimension)
    {
        switch (dimension)
        {
            case 0:
                return this->d1;
            case 1:
                return this->d2;
            case 2:
                return this->d3;
            default:
                printf("Error in List3d getSize()");
                throw;
        }
    }
private:
    size_t d1, d2, d3;
    std::vector<T> data;
};