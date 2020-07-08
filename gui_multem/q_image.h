#ifndef Q_IMAGE_H
#define Q_IMAGE_H

#include <thrust/complex.h>
#include "types.cuh"

#include <QWidget>
#include <QtWidgets>
#include "q_colormap.h"

class Data_Image : public QWidget
{
	Q_OBJECT

public:
	using TVector_s = thrust::host_vector<float>;

	explicit Data_Image(QWidget *parent = Q_NULLPTR)
	: QWidget(parent), nx(0), ny(0), nxy(0), image(nullptr),
	im_zoom_0(1.0), im_zoom(1.0), im_zoom_d(1.25)
	{
		im_zoom_min = 0.1*im_zoom_0;
		im_zoom_max = 20.0*im_zoom_0;

		im_pos = QPointF(0, 0);
		im_pos_press = im_pos;
		window_size = QSize(0, 0);

		FileName = "";

		setContextMenuPolicy(Qt::CustomContextMenu);

		mn_Contex = new QMenu();
		act_save_as_image = mn_Contex->addAction("Save as image");
		act_save_as_binary = mn_Contex->addAction("Save as Binary");

		connect(this, SIGNAL(customContextMenuRequested(const QPoint &)),SLOT(custom_context_menu(const QPoint &)));
		connect(act_save_as_image, SIGNAL(triggered()),SLOT(save_as_tif()));
		connect(act_save_as_binary, SIGNAL(triggered()),SLOT(save_as_binary()));
	}
	
	~Data_Image(){}

	QSize minimumSizeHint() const Q_DECL_OVERRIDE
	{
		return sizeHint()/4;
	}

	QSize sizeHint() const Q_DECL_OVERRIDE
	{
		auto size_xy = qRound(max(im_zoom_0*nx, im_zoom_0*ny));
		return QSize(size_xy, size_xy);
	}

	QSize data_size()
	{
		return QSize(nx, ny);
	}

	void set_data_size(int nx_i, int ny_i)
	{
		nx = nx_i;
		ny = ny_i;
		nxy = nx*ny;

		delete image;
		image = new QImage(nx, ny, QImage::Format_RGB32);

		set_default_zoom_pos();
	}

	void set_window_size(QSize ws)
	{
		window_size = ws;
		set_default_zoom_pos();
	}

	void set_default_zoom_pos()
	{
		auto data_size_f = QPointF(data_size().width(), data_size().height());
		im_zoom_0 = qreal(window_size.height())/qreal(max(data_size_f.rx(), data_size_f.ry()));
		im_zoom = im_zoom_0;
		im_zoom_min = 0.1*im_zoom_0;
		im_zoom_max = 20.0*im_zoom_0;

		im_pos.setX(0.5*(window_size.width() - im_zoom_0*data_size_f.rx()));
		im_pos.setY(0.5*(window_size.height() - im_zoom_0*data_size_f.ry()));
	}

	void scale(qreal factor)
	{
		im_zoom = factor;
		update();
	}

	void pos_scale(QPointF pos, qreal factor)
	{
		im_pos = pos;
		im_zoom = factor;
		update();
	}

	void draw(TVector_s *data, eColormap colormap, eScale_Type scale_type, double scale_factor)
	{
		Colormap map(colormap);

		TVector_s &data_s = *data;
		auto data_s_minmax = std::minmax_element(data_s.begin(), data_s.end());
		double data_s_min = *(data_s_minmax.first);
		double data_s_max = *(data_s_minmax.second);

		switch(scale_type)
		{
			case eSLT_Linear:
			{
				double vs_min = data_s_min;
				double vs_max = data_s_max;
				double y_m = 255/(vs_max-vs_min);

				for(auto iy = 0; iy<ny; iy++)
				{
					auto line = (QRgb*)(image->scanLine(iy));
					for(auto ix = 0; ix<nx; ix++)
					{
						auto v_r = data_s[ix*ny+iy];
						auto v_u = real_2_uchar(v_r, vs_min, y_m);
						line[ix] = map(v_u);
					}
				}
			}
			break;
			case eSLT_Log:
			{
				double vs_min = log_scale(data_s_min);
				double vs_max = log_scale(data_s_max);
				double y_m = 255/(vs_max-vs_min);

				for(auto iy = 0; iy<ny; iy++)
				{
					auto line = (QRgb*)(image->scanLine(iy));
					for(auto ix = 0; ix<nx; ix++)
					{
						auto v_r = log_scale(data_s[ix*ny+iy]);
						auto v_u = real_2_uchar(v_r, vs_min, y_m);
						line[ix] = map(v_u);
					}
				}
			}
			break;
			case eSLT_Power:
			{
				double vs_min = power_scale(data_s_min, scale_factor);
				double vs_max = power_scale(data_s_max, scale_factor);
				double y_m = 255/(vs_max-vs_min);

				for(auto iy = 0; iy<ny; iy++)
				{
					auto line = (QRgb*)(image->scanLine(iy));
					for(auto ix = 0; ix<nx; ix++)
					{
						auto v_r = power_scale(data_s[ix*ny+iy], scale_factor);
						auto v_u = real_2_uchar(v_r, vs_min, y_m);
						line[ix] = map(v_u);
					}
				}
			}
			break;
		}

		update();
	}

signals:
	void sg_save_as_binary();

public slots:
	void set_normalSize()
	{
		set_default_zoom_pos();
		update();
	}

	void custom_context_menu(const QPoint &point)
	{
		mn_Contex->popup(mapToGlobal(point));
	}

	void save_data(QString name)
	{
		image->save(name);
	}

	void save_as_tif()
	{
		FileName = QFileDialog::getSaveFileName(0,"Save file", QDir::currentPath(),
		"PNG - Portable Network Graphics (*.png);TIF -  TIFF Bitmap (*.tif)",
		new QString("PNG - Portable Network Graphics (*.png)"));

		if (!FileName.isEmpty())
		{
		   image->save(FileName);
		}
	}

	void save_as_binary()
	{
		emit sg_save_as_binary();
	}

protected:
	void paintEvent(QPaintEvent * /*event*/) Q_DECL_OVERRIDE
	{
		QPainter painter(this);
		painter.setRenderHint(QPainter::Antialiasing, false);
		painter.setRenderHint(QPainter::SmoothPixmapTransform, false);
		painter.setRenderHint(QPainter::HighQualityAntialiasing, false);
		painter.translate(im_pos);
		painter.scale(im_zoom, im_zoom);
		painter.drawImage(QPointF(0, 0), *image);
	}
	
	void wheelEvent(QWheelEvent * event) Q_DECL_OVERRIDE
	{
		qreal delta = (event->delta()<0)?1.0/im_zoom_d:im_zoom_d;
		auto zoom_t = im_zoom*delta;

		if((im_zoom_min<zoom_t)&&(zoom_t<im_zoom_max))
		{
			im_pos = delta*im_pos + (1.0-delta)*event->posF();
			im_zoom *= delta;
		}

		update();
	}
	
	void mouseMoveEvent(QMouseEvent *event) Q_DECL_OVERRIDE
	{
		event->accept();
//		if(event->buttons()==Qt::LeftButton)
//		{
//			im_pos = im_pos + (event->pos() - im_pos_press)*im_zoom;
//			update();
//		}
		emit mouseMove(event->pos());
	}
	
	void mousePressEvent(QMouseEvent *event) Q_DECL_OVERRIDE
	{
		event->accept();
		im_pos_press = event->pos();
		emit mousePress(event->pos());
	}
	
	void mouseReleaseEvent(QMouseEvent *event) Q_DECL_OVERRIDE
	{
		event->accept();
		im_pos_press = QPointF(0, 0);
		emit mouseRelease(event->pos());
	}

	void mouseDoubleClickEvent(QMouseEvent *event) Q_DECL_OVERRIDE
	{
		event->accept();
		set_normalSize();
		emit mouseDoubleClick();
	}

signals:
	void mouseDoubleClick();
	void mouseMove(QPoint Point);
	void mousePress(QPoint Point);
	void mouseRelease(QPoint Point);

private:
	QString FileName;

	double log_scale(double x)
	{
		const double ee = 1e-6;
		return log(ee+x);
	}

	double power_scale(double x, double factor)
	{
		return pow(x, factor);
	}

	uchar real_2_uchar(double x, double x_min, double y_m)
	{
		auto v = (uchar)(y_m*(x-x_min));
		return qBound<uchar>(0, v, 255);
	}

	int nx;
	int ny;
	int nxy;

	QImage *image;
	QSize window_size;

	qreal im_zoom_0;		// Initial zoom factor
	qreal im_zoom;		  // zoom factor
	qreal im_zoom_min;	  // zoom factor
	qreal im_zoom_max;	  // zoom factor
	qreal im_zoom_d;		// step
	QPointF im_pos;		 //
	QPointF im_pos_press;   //

	QMenu *mn_Contex;
	QAction *act_save_as_image;
	QAction *act_save_as_binary;
};

#endif // Q_IMAGE_H
