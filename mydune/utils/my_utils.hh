#ifndef MY_UTILS_HH
#define MY_UTILS_HH

#include <vector>

class QPoint
{
    public:
        QPoint();
        QPoint(double x0,double y0,double z0, double q0)
        {
            x=x0;
            y=y0;
            z=z0;
            q=q0;
        }
        ~QPoint() {}
        double x,y,z,q;
};


class QCluster : public std::vector<QPoint>
{
public:
    int sliceID = -1;

    QCluster() = default;

    QCluster(const QCluster& other)
        : std::vector<QPoint>(other), sliceID(other.sliceID)
    {}

    QCluster& operator=(const QCluster& other)
    {
        if (this != &other)
        {
            std::vector<QPoint>::operator=(other);
            sliceID = other.sliceID;
        }
        return *this;
    }
};




#endif