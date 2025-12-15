#ifndef MY_UTILS_HH
#define MY_UTILS_HH

#include <vector>

class QPoint
{
    public:
        QPoint()
        {
        }
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
        int objID = -1;

        QCluster() = default;

        QCluster(const QCluster& other)
            : std::vector<QPoint>(other), objID(other.objID)
        {}

        QCluster& operator=(const QCluster& other)
        {
            if (this != &other)
            {
                std::vector<QPoint>::operator=(other);
                objID = other.objID;
            }
            return *this;
        }
};


class QFlash
{
    public: 
        QFlash(){};
        ~QFlash(){};
        int flashID = -1;
        std::vector<double> PE_CH;
        double y,z,y_err,z_err;
        double time, time_err;
};




#endif