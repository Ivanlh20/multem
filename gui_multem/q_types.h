#ifndef QT_TYPES_H
#define QT_TYPES_H

#include <limits>
#include <QColor>
#include <QString>
#include <QByteArray>

class QLineEdit;
class QComboBox;
class QWidget;

/************************Input File*******************************/
enum eInput_File
{
	eIF_txt = 1, eIF_pdb = 2
};

/************************Color map*******************************/
enum eColormap
{
	eCM_Gray = 1, eCM_Hot = 2, eCM_Cool = 3, eCM_Jet = 4,
	eCM_Copper = 5, eCM_Summer = 6, eCM_Autumn = 7, eCM_Winter = 8
};

/**************************Scale*********************************/
enum eScale_Type
{
	eSLT_Linear = 1, eSLT_Power = 2, eSLT_Log = 3
};


const QString Ags = QString::fromWCharArray(L"\x00C5");
const QString rAgs = QString::fromWCharArray(L"\x00C5\x207B\x00B9");
const QString elec = QString::fromWCharArray(L"e\x207B");
const double double_lim = std::numeric_limits<double>::max()/2;
const double int_lim = std::numeric_limits<int>::max()/2;

inline
QString sups_mone()
{
	return QString::fromWCharArray(L"\x207B\x00B9");
}

inline
QByteArray c_nm_to_qba(int nm, const char *units)
{
	QString str = "C<sub style=""font-family:calibri;font-size:x-large;"">";
	str += QString::number(nm) + "</sub>" + QString::fromUtf8(units);
	return str.toUtf8();
}

inline
QByteArray phi_nm_to_qba(int nm, const char *units)
{
	QString str = "<span style=""font-family:symbol;font-size:large"">f</span><sub style=""font-family:calibri;font-size:x-large;"">";
	str += QString::number(nm) + "</sub>" + QString::fromUtf8(units);
	return str.toUtf8();
}

inline
QByteArray x_sub_to_qba(const char *x, const char *s, QString fs="large")
{
	QString str = QString::fromUtf8(x) + "<sub style=""font-size:" + fs + """>";
	str += QString::fromUtf8(s) + "</sub>";
	return str.toUtf8();
}

inline
QByteArray n_atoms_types_qba(int na, int nt)
{
	QString str_natoms = "N<sub style=""font-size:large"">atoms</sub> = " + QString::number(na);
	QString str_ntypes = "N<sub style=""font-size:large"">atom types</sub> = " + QString::number(nt);
	QString str = str_natoms + "   -   " + str_ntypes;
	return str.toUtf8();
}

inline
QByteArray n_slices_qba(int na)
{
	QString str_nslices = "N<sub style=""font-size:large"">slices</sub> = " + QString::number(na);
	return str_nslices.toUtf8();
}

inline
QString thk_value_qstr(double value)
{
	QString str_thk = " " + QString::number(value, 'f', 4) + " " + Ags + " ";
	return str_thk;
}

inline
QString det_value_qstr(double lambda, double g_inner, double g_outer)
{
	g_inner = 1000*asin(lambda*g_inner);
	g_outer = 1000*asin(lambda*g_outer);
	QString str_det = " [" + QString::number(g_inner, 'f', 3) + " - " + QString::number(g_outer, 'f', 3) + "] mrad ";
	return str_det;
}

template<class T>
QVariant enum_to_QVar(const T &data)
{
	return QVariant(static_cast<int>(data));
}

template<class T>
T QVar_to_enum(const QVariant &data)
{
	return static_cast<T>(data.toInt());
}

inline
void disable_item_QComboBox(QComboBox *cb, int index, bool disabled=true)
{
	QStandardItemModel *model = qobject_cast<QStandardItemModel*>(cb->model());
	QStandardItem *item = model->item(index);
	item->setFlags(disabled? item->flags() & ~Qt::ItemIsEnabled:
							 item->flags() | Qt::ItemIsEnabled);
}

inline
void set_non_editable_QLineEdit(QLineEdit *le)
{
	le->setReadOnly(true);
	le->setEnabled(false);

//	auto palette = new QPalette();
//	palette->setColor(QPalette::Base,Qt::gray);
//	palette->setColor(QPalette::Text,Qt::darkGray);

//	le->setReadOnly(true);
//	le->setPalette(*palette);
}

#endif // QT_TYPES_H
