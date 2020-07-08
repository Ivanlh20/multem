#ifndef Colormap_H
#define Colormap_H

#include <vector>
#include <algorithm>
#include <QColor>
#include <QVector3D>
#include "q_types.h"

class Colormap
{
public:
	Colormap(const eColormap &colormap)
	{
		points.resize(128);
		map.resize(256);

		set_colormap(colormap);
	}

	Colormap()
	{
		points.resize(128);
		map.resize(256);

		clear();
	}

	~Colormap(){}

	void set_colormap(const eColormap &colormap)
	{
		const int m = 256;
		int ri, gi, bi;

		clear();

		switch(colormap)
		{
			case eCM_Gray:
			{
				for(auto ik=0; ik<m; ik++)
				{
					double r = double(ik)/double(m-1);
					ri = gi = bi = qBound(0, qRound((m-1)*r), m-1);
					map[ik] = qRgb(ri, gi, bi);
				}
			}
			break;
			case eCM_Cool:
			{
				for(auto ik=0; ik<m; ik++)
				{
					double r = double(ik)/double(m-1);
					ri = qBound(0, qRound((m-1)*r), m-1);
					gi = (m-1) - ri;
					bi = m-1;
					map[ik] = qRgb(ri, gi, bi);
				}
			}
			break;
			case eCM_Hot:
			{
				auto n = qRound(3.0*m/8.0);
				for(auto ik=0; ik<m; ik++)
				{
					double r = (ik<=n)?ik/double(n):1.0;
					double g = (ik<=n)?0.0:(ik<=2*n)?double(ik-n)/double(n):1.0;
					double b = (ik<=2*n)?0.0:double(ik-2*n)/double(m-2*n);

					ri = qBound(0, qRound((m-1)*r), m-1);
					gi = qBound(0, qRound((m-1)*g), m-1);
					bi = qBound(0, qRound((m-1)*b), m-1);

					map[ik] = qRgb(ri, gi, bi);
				}
			}
			 break;
			case eCM_Jet:
			{
				auto n = qRound(m/4.0);
				for(auto ik=0; ik<m; ik++)
				{
					ri = ik-1.5*n;
					double r = (ri<=n)?double(ri)/double(n):(ri<=2*n)?1.0:(ri<=3*n)?double(3*n-ri)/double(n):0.0;
					gi = ik-0.5*n;
					double g = (gi<=n)?double(gi)/double(n):(gi<=2*n)?1.0:(gi<=3*n)?double(3*n-gi)/double(n):0.0;
					bi = ik+0.5*n;
					double b = (bi<=n)?double(bi)/double(n):(bi<=2*n)?1.0:(bi<=3*n)?double(3*n-bi)/double(n):0.0;
					ri = qBound(0, qRound((m-1)*r), m-1);
					gi = qBound(0, qRound((m-1)*g), m-1);
					bi = qBound(0, qRound((m-1)*b), m-1);

					map[ik] = qRgb(ri, gi, bi);
				}
			}
			break;
			case eCM_Copper:
			{
				for(auto ik=0; ik<m; ik++)
				{
					double t = double(ik)/double(m-1);
					double r = qMin(1.0, 1.250*t);
					double g = qMin(1.0, 0.7812*t);
					double b = qMin(1.0, 0.4975*t);

					ri = qBound(0, qRound((m-1)*r), m-1);
					gi = qBound(0, qRound((m-1)*g), m-1);
					bi = qBound(0, qRound((m-1)*b), m-1);

					map[ik] = qRgb(ri, gi, bi);
				}
			}
			break;
			case eCM_Summer:
			{
				for(auto ik=0; ik<m; ik++)
				{
					double r = double(ik)/double(m-1);
					double g = 0.5*(1.0+r);
					double b = 0.4;

					ri = qBound(0, qRound((m-1)*r), m-1);
					gi = qBound(0, qRound((m-1)*g), m-1);
					bi = qBound(0, qRound((m-1)*b), m-1);

					map[ik] = qRgb(ri, gi, bi);
				}
			}
			break;
			case eCM_Autumn:
			{
				for(auto ik=0; ik<m; ik++)
				{
					double r = 1.0;
					double g = double(ik)/double(m-1);
					double b = 0.0;

					ri = qBound(0, qRound((m-1)*r), m-1);
					gi = qBound(0, qRound((m-1)*g), m-1);
					bi = qBound(0, qRound((m-1)*b), m-1);

					map[ik] = qRgb(ri, gi, bi);
				}
			}
			break;
			case eCM_Winter:
			{
				for(auto ik=0; ik<m; ik++)
				{
					double r = 0.0;
					double g = double(ik)/double(m-1);
					double b = 1.0-0.5*r;

					ri = qBound(0, qRound((m-1)*r), m-1);
					gi = qBound(0, qRound((m-1)*g), m-1);
					bi = qBound(0, qRound((m-1)*b), m-1);

					map[ik] = qRgb(ri, gi, bi);
				}
			}
			break;
		}
	}

	QRgb operator()(const uchar &v)
	{
		return map[v];
	}

private:
	int npoints;
	std::vector<QRgb> points;
	std::vector<QRgb> map;

	void clear()
	{
		std::fill(points.begin(), points.end(), qRgba(0, 0, 0, 0));
		std::fill(map.begin(), map.end(), qRgba(0, 0, 0, 0));
	}

	void getLinearmap()
	{
		if(points.empty()) return;

		int ipoints;
		int r, g, b;
		QVector3D V0, V1, u0;
		int idx0, idx1, idx;

		// get initial point
		idx0 = qAlpha(points[0]);
		V0.setX(qRed(points[0]));
		V0.setY(qGreen(points[0]));
		V0.setZ(qBlue(points[0]));

		for(ipoints=1; ipoints<points.size(); ipoints++){
			idx1 = qAlpha(points[ipoints]);
			V1.setX(qRed(points[ipoints]));
			V1.setY(qGreen(points[ipoints]));
			V1.setZ(qBlue(points[ipoints]));
			u0 = (V1-V0)/(idx1-idx0);

			for(idx=idx0; idx<idx1; idx++){
				V1 = V0 + u0*(idx-idx0);
				r = qBound(0, qRound(V1.x()), 255);
				g = qBound(0, qRound(V1.y()), 255);
				b = qBound(0, qRound(V1.z()), 255);
				map[idx] = qRgb(r, g, b);
			}

			V0 =  V1;
		}
		map[idx1] = qRgb(qRed(points.back()), qGreen(points.back()), qBlue(points.back()));
	}
};

#endif // Colormap_H
