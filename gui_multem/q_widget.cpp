#include "q_widget.h"

#include <QAction>
#include <QtEvents>
#include <QFrame>
#include <QMainWindow>
#include <QPainter>
#include <QColor>
#include <QSpinBox>
#include <QTabWidget>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QLineEdit>
#include <QComboBox>
#include <QParallelAnimationGroup>
#include <QScrollArea>
#include <QToolButton>
#include <QPropertyAnimation>

QColor from_string_to_color(const QString &name)
{
	if (name == "Black")
		return QColor("#D8D8D8");
	if (name == "White")
		return QColor("#F1F1F1");
	if (name == "Red")
		return QColor("#F1D8D8");
	if (name == "Green")
		return QColor("#D8E4D8");
	if (name == "Blue")
		return QColor("#D8D8F1");
	if (name == "Yellow")
		return QColor("#F1F0D8");

	return QColor(name).light(110);
}

/*************************** MColorMap ****************************/

void MColorMap::set(Spec color_map_i)
{
	colormap = color_map_i;

	int v_max = qint32(map.size())-1;

	std::fill(map.begin(), map.end(), qRgba(0, 0, 0, 0));

	auto to_int = [&v_max](const qreal &v_n)
	{
		return qBound<int>(0, qRound(v_max*v_n), v_max);
	};

	switch(colormap)
	{
	case eCM_Gray:
	{
		for(auto iv=0; iv<map.size(); iv++)
		{
			auto r = qreal(iv)/qreal(v_max);

			auto r_i = to_int(r);
			auto g_i = r_i;
			auto b_i = r_i;

			map[iv] = qRgb(r_i, g_i, b_i);
		}
	}
		break;
	case eCM_Cool:
	{
		for(auto iv=0; iv<map.size(); iv++)
		{
			auto r = qreal(iv)/qreal(v_max);

			auto r_i = to_int(r);
			auto g_i = v_max - r_i;
			auto b_i = v_max;

			map[iv] = qRgb(r_i, g_i, b_i);
		}
	}
		break;
	case eCM_Hot:
	{
		auto v_lim = qRound(3*(v_max+1)/8.0);

		for(auto iv=0; iv<map.size(); iv++)
		{
			auto r = (iv<=v_lim)?(qreal(iv)/qreal(v_lim)):1.0;
			auto g = (iv<=v_lim)?0.0:(iv<=2*v_lim)?(qreal(iv-v_lim)/qreal(v_lim)):1.0;
			auto b = (iv<=2*v_lim)?0.0:(qreal(iv-2*v_lim)/qreal(v_max-2*v_lim));

			auto r_i = to_int(r);
			auto g_i = to_int(g);
			auto b_i = to_int(b);

			map[iv] = qRgb(r_i, g_i, b_i);
		}
	}
		break;
	case eCM_Jet:
	{
		auto v_lim = qRound((v_max+1)/4.0);

		for(auto iv=0; iv<map.size(); iv++)
		{
			auto r = iv-1.5*qreal(v_lim);
			r = (r<=v_lim)?qreal(r)/qreal(v_lim):(r<=2*v_lim)?1.0:(r<=3*v_lim)?qreal(3*v_lim-r)/qreal(v_lim):0.0;
			auto g = iv-0.5*v_lim;
			g = (g<=v_lim)?qreal(g)/qreal(v_lim):(g<=2*v_lim)?1.0:(g<=3*v_lim)?qreal(3*v_lim-g)/qreal(v_lim):0.0;
			auto b = iv+0.5*v_lim;
			b = (b<=v_lim)?qreal(b)/qreal(v_lim):(b<=2*v_lim)?1.0:(b<=3*v_lim)?qreal(3*v_lim-b)/qreal(v_lim):0.0;

			auto r_i = to_int(r);
			auto g_i = to_int(g);
			auto b_i = to_int(b);

			map[iv] = qRgb(r_i, g_i, b_i);
		}
	}
		break;
	case eCM_Copper:
	{
		for(auto iv=0; iv<map.size(); iv++)
		{
			auto t = qreal(iv)/qreal(v_max);
			auto r = qMin<qreal>(1.0, 1.250*t);
			auto g = qMin<qreal>(1.0, 0.7812*t);
			auto b = qMin<qreal>(1.0, 0.4975*t);

			auto r_i = to_int(r);
			auto g_i = to_int(g);
			auto b_i = to_int(b);

			map[iv] = qRgb(r_i, g_i, b_i);
		}
	}
		break;
	case eCM_Summer:
	{
		for(auto iv=0; iv<map.size(); iv++)
		{
			auto r = qreal(iv)/qreal(v_max);
			auto g = 0.5*(1.0+r);
			auto b = 0.4;

			auto r_i = to_int(r);
			auto g_i = to_int(g);
			auto b_i = to_int(b);

			map[iv] = qRgb(r_i, g_i, b_i);
		}
	}
		break;
	case eCM_Autumn:
	{
		for(auto iv=0; iv<map.size(); iv++)
		{
			auto r = 1.0;
			auto g = qreal(iv)/qreal(v_max);
			auto b = 0.0;

			auto r_i = to_int(r);
			auto g_i = to_int(g);
			auto b_i = to_int(b);

			map[iv] = qRgb(r_i, g_i, b_i);
		}
	}
		break;
	case eCM_Winter:
	{
		for(auto iv=0; iv<map.size(); iv++)
		{
			auto r = 0.0;
			auto g = qreal(iv)/qreal(v_max);
			auto b = 1.0-0.5*r;

			auto r_i = to_int(r);
			auto g_i = to_int(g);
			auto b_i = to_int(b);

			map[iv] = qRgb(r_i, g_i, b_i);
		}
	}
		break;
	}
}

/***************************** MGroupBox_BE *****************************/

MGroupBox_BE::MGroupBox_BE(const QString &title, QWidget *parent)
	: QWidget(parent)
{
	header_title = new QGroupBox(title);
	header_title->setAlignment(Qt::AlignHCenter);
	header_title->setFlat(true);
	header_title->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Maximum);

	expand_buttom = new QToolButton;
	expand_buttom->setStyleSheet("QToolButton { border: none; }");
	expand_buttom->setToolButtonStyle(Qt::ToolButtonIconOnly);
	expand_buttom->setArrowType(Qt::ArrowType::LeftArrow);
	expand_buttom->setCheckable(true);
	expand_buttom->setChecked(false);

	/********************************************************/
	content_frame = new QFrame;

	main_layout = new QGridLayout;
	main_layout->setVerticalSpacing(0);
	main_layout->setContentsMargins(0, 0, 0, 0);
	main_layout->addWidget(header_title, 0, 0);
	main_layout->addWidget(expand_buttom, 0, 1, 1, 1, Qt::AlignRight);
	main_layout->addWidget(content_frame, 1, 0, 1, 2);

	setLayout(main_layout);

	QObject::connect(expand_buttom, &QToolButton::toggled, [this](const bool checked) {
		expand_buttom->setArrowType(checked ? Qt::ArrowType::DownArrow : Qt::ArrowType::LeftArrow);
		content_frame->setVisible(checked);
	});
	content_frame->setVisible(expand_buttom->isChecked());
}
