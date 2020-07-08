#ifndef MWIDGET_H
#define MWIDGET_H

#include <vector>

#include <QWidget>
#include <QFrame>
#include <QDockWidget>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QLabel>
#include <QGroupBox>
#include <QToolButton>

class MColorMap
{
public:
	enum Spec
	{
		eCM_Invalid = 0, eCM_Gray = 1, eCM_Hot = 2, eCM_Cool = 3, eCM_Jet = 4,
		eCM_Copper = 5, eCM_Summer = 6, eCM_Autumn = 7, eCM_Winter = 8
	};

	inline
	MColorMap(Spec color_map_i=eCM_Gray):map(256, qRgba(0, 0, 0, 0))
	{
		set(color_map_i);
	}

	inline
	QRgb operator()(const uchar &iv) const
	{
		return map[iv];
	}

	void set(Spec color_map_i);

	inline
	Spec get() const
	{
		return colormap;
	}

private:
	std::vector<QRgb> map;
	Spec colormap;
};

/***************************** MLine_H *****************************/
class MLine_H: public QFrame
{
	Q_OBJECT
public:
	explicit MLine_H(QWidget *parent = Q_NULLPTR): QFrame(parent)
	{
		setFrameStyle(QFrame::HLine | QFrame::Sunken);
		setLineWidth(1);
	}
};

/***************************** MFrame_V *****************************/
class MFrame_V: public QFrame
{
	Q_OBJECT
public:
	explicit MFrame_V(QWidget *parent = Q_NULLPTR)
		: QFrame(parent), layout(new QVBoxLayout),
		  szHint(-1, -1), minSzHint(50, 50)
	{
		setFrameStyle(QFrame::Box | QFrame::Sunken);
		setLineWidth(1);
		setMidLineWidth(1);
		layout->setContentsMargins(1, 2, 1, 0);
		//layout->setVerticalSpacing(2);
		layout->setSizeConstraint(QLayout::SetMinimumSize);
		setLayout(layout);
	}

	QSize sizeHint() const Q_DECL_OVERRIDE
	{
		return szHint;
	}

	QSize minimumSizeHint() const Q_DECL_OVERRIDE
	{
		return minSzHint;
	}

	QVBoxLayout *layout;
private:
	QSize szHint;
	QSize minSzHint;
};

/***************************** MFrame_G *****************************/

class MFrame_G: public QFrame
{
	Q_OBJECT
public:
	explicit MFrame_G(QWidget *parent = Q_NULLPTR)
		: QFrame(parent), layout(new QGridLayout)
		, szHint(-1, -1), minSzHint(50, 50)
	{
		setFrameStyle(QFrame::Box | QFrame::Sunken);
		setLineWidth(1);
		setMidLineWidth(1);

		layout->setContentsMargins(1, 2, 1, 0);
		layout->setVerticalSpacing(2);
		layout->setSizeConstraint(QLayout::SetMinimumSize);
		setLayout(layout);
	}

	QSize sizeHint() const Q_DECL_OVERRIDE
	{
		return szHint;
	}

	QSize minimumSizeHint() const Q_DECL_OVERRIDE
	{
		return minSzHint;
	}

	QGridLayout *layout;
private:
	QSize szHint;
	QSize minSzHint;
};

/***************************** MGroupBox_BE *****************************/
class MGroupBox_BE : public QWidget
{
	Q_OBJECT

public:
	explicit MGroupBox_BE(const QString & title = "", QWidget *parent = Q_NULLPTR);

	void setExpand(bool expand)
	{
		expand_buttom->setChecked(expand);
	}

protected:
	QGridLayout *main_layout;
	QGroupBox *header_title;
	QToolButton *expand_buttom;
	QFrame *content_frame;
};

/***************************** MGroupBox_HE *****************************/

class MGroupBox_HE : public MGroupBox_BE
{
	Q_OBJECT

public:
	explicit MGroupBox_HE(const QString & title = "", QWidget *parent = Q_NULLPTR)
		:MGroupBox_BE::MGroupBox_BE(title, parent), layout(new QHBoxLayout)
	{
		layout->setContentsMargins(2, 1, 1, 0);
		this->content_frame->setLayout(layout);
	}

	void addWidget(QWidget *widget)
	{
		layout->addWidget(widget);
	}

	template<class... TArgs>
	void addWidget(QWidget *widget, TArgs... args)
	{
		layout->addWidget(widget);
		addWidget(args...);
	}

	QHBoxLayout *layout;
};

/***************************** MGroupBox_VE *****************************/

class MGroupBox_VE : public MGroupBox_BE
{
	Q_OBJECT

public:
	explicit MGroupBox_VE(const QString & title = "", QWidget *parent = Q_NULLPTR)
		:MGroupBox_BE::MGroupBox_BE(title, parent), layout(new QVBoxLayout)
	{
		layout->setSpacing(0);
		layout->setContentsMargins(2, 1, 1, 0);
		content_frame->setLayout(layout);
	}

	void addWidget(QWidget *widget)
	{
		layout->addWidget(widget);
	}

	template<class... TArgs>
	void addWidget(QWidget *widget, TArgs... args)
	{
		layout->addWidget(widget);
		addWidget(args...);
	}

	QVBoxLayout *layout;
};

/***************************** MGroupBox_GE *****************************/

class MGroupBox_GE : public MGroupBox_BE
{
	Q_OBJECT

public:
	explicit MGroupBox_GE(const QString & title = "", QWidget *parent = Q_NULLPTR)
		:MGroupBox_BE::MGroupBox_BE(title, parent), layout(new QGridLayout)
	{
		layout->setVerticalSpacing(1);
		layout->setHorizontalSpacing(2);
		layout->setContentsMargins(2, 1, 1, 0);
		this->content_frame->setLayout(layout);
	}

	void addWidget(QWidget *widget, int row, int column, Qt::Alignment alignment = Qt::Alignment())
	{
		layout->addWidget(widget, row, column, alignment);
	}

	void addWidget(QWidget *widget, int fromRow, int fromColumn, int rowSpan, int columnSpan, Qt::Alignment alignment = Qt::Alignment())
	{
		layout->addWidget(widget, fromRow, fromColumn, rowSpan, columnSpan, alignment);
	}

	QGridLayout *layout;
};


/***************************** MGroupBox_H *****************************/

class MGroupBox_H: public QGroupBox
{
	Q_OBJECT
public:
	explicit MGroupBox_H(const QString &title = "", QWidget *parent = Q_NULLPTR)
		: QGroupBox(title, parent), layout(new QHBoxLayout)
	{
		setAlignment(Qt::AlignHCenter);
		setFlat(true);
		layout->setSpacing(1);
		layout->setContentsMargins(1, 2, 1, 2);
		setLayout(layout);
	}

	void addWidget(QWidget *widget)
	{
		layout->addWidget(widget);
	}

	template<class... TArgs>
	void addWidget(QWidget *widget, TArgs... args)
	{
		layout->addWidget(widget);
		addWidget(args...);
	}

	QHBoxLayout *layout;
};

/***************************** MGroupBox_V *****************************/

class MGroupBox_V: public QGroupBox
{
	Q_OBJECT
public:
	explicit MGroupBox_V(const QString &title = "", QWidget *parent = Q_NULLPTR)
		: QGroupBox(title, parent), layout(new QVBoxLayout)
	{
		setAlignment(Qt::AlignHCenter);
		setFlat(true);
		layout->setSpacing(1);
		layout->setContentsMargins(1, 2, 1, 2);
		setLayout(layout);
	}

	void addWidget(QWidget *widget)
	{
		layout->addWidget(widget);
	}

	template<class... TArgs>
	void addWidget(QWidget *widget, TArgs... args)
	{
		layout->addWidget(widget);
		addWidget(args...);
	}

	QVBoxLayout *layout;
};

/***************************** MGroupBox_G *****************************/

class MGroupBox_G: public QGroupBox
{
	Q_OBJECT
public:
	explicit MGroupBox_G(const QString &title = "", QWidget *parent = Q_NULLPTR)
		: QGroupBox(title, parent), layout(new QGridLayout)
	{
		setAlignment(Qt::AlignHCenter);
		setFlat(true);
		layout->setVerticalSpacing(0);
		layout->setHorizontalSpacing(2);
		layout->setContentsMargins(1, 2, 1, 2);
		setLayout(layout);
	}

	void addWidget(QWidget *widget, int row, int column, Qt::Alignment alignment = Qt::Alignment())
	{
		layout->addWidget(widget, row, column, alignment);
	}

	void addWidget(QWidget *widget, int fromRow, int fromColumn, int rowSpan, int columnSpan, Qt::Alignment alignment = Qt::Alignment())
	{
		layout->addWidget(widget, fromRow, fromColumn, rowSpan, columnSpan, alignment);
	}

	QGridLayout *layout;
};

#endif // MWIDGET_H
