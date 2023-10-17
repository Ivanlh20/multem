#ifndef DATA_VIEWER_H
#define DATA_VIEWER_H

#include "types.cuh"
#include "output_multislice.hpp"

#include <QWidget>
#include <QtWidgets>

#include "q_types.h"
#include "q_image.h"

class Data_Viewer: public QDialog
{
	Q_OBJECT
public:
	template<class TOutput_Multislice>
	explicit Data_Viewer(TOutput_Multislice &output_multislice_i, QDialog *parent = Q_NULLPTR): QDialog(parent)
	{
		setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
		MouSetPos = QPoint(0,0);

		QDesktopWidget *desktopWidget = QApplication::desktop();
		const QRect screenGeometry = desktopWidget->screenGeometry(desktopWidget->primaryScreen());
		window_size = screenGeometry.size()*3/5;
		auto side_min = qMin(window_size.width(), window_size.height());
		window_size.setWidth(side_min);
		window_size.setHeight(side_min);

		/******************************************************/
		create_wigets();

		/******************************************************/
		output_multislice = output_multislice_i;

		image->set_data_size(output_multislice.nx, output_multislice.ny);
		image->set_window_size(window_size);

		/******************************************************/
		sl_thickness->setMinimum(1);
		sl_thickness->setMaximum(output_multislice.thick.size());
		sl_thickness->setValue(1);

		if(is_thickness_available())
		{
		   lb_thickness_value->setText(thk_value_qstr(output_multislice.thick[0]));
		}

		/******************************************************/
		sl_detector->setMinimum(1);
		sl_detector->setMaximum(output_multislice.ndetector);
		sl_detector->setValue(1);

		if(output_multislice.is_STEM())
		{
			auto lambda = mt::get_lambda(output_multislice.E_0);
			auto g_inner = output_multislice.detector.g_inner[0];
			auto g_outer = output_multislice.detector.g_outer[0];

			lb_detector_value->setText(det_value_qstr(lambda, g_inner, g_outer));
		}

		/******************************************************/
		activate_options();

		/******************************************************/
		connect(image, SIGNAL(mouseDoubleClick()), this, SLOT(set_normalSize()));

//		connect(image_tot, SIGNAL(mouseMove(QPoint)), this, SLOT(qt_image_mouseMoveEvent(QPoint)));
//		connect(image_tot, SIGNAL(mousePress(QPoint)), this, SLOT(qt_image_mousePressEvent(QPoint)));
//		connect(image_tot, SIGNAL(mouseReleaseEvent(QPoint)), this, SLOT(qt_image_mouseReleaseEvent(QPoint)));
//		connect(cb_colormap, SIGNAL(currentIndexChanged(int)), this, SLOT(on_cb_ColorMap_currentIndexChanged(int)));


		connect(cb_show_ctype, SIGNAL(currentIndexChanged(int)), SLOT(cb_show_ctype_currentIndexChanged(int)));
		connect(cb_colormap, SIGNAL(currentIndexChanged(int)), SLOT(cb_colormap_currentIndexChanged(int)));

		connect(cb_scale, SIGNAL(currentIndexChanged(int)), SLOT(cb_scale_currentIndexChanged(int)));

		connect(ckb_block_g0, SIGNAL(stateChanged(int)), SLOT(ckb_block_g0_stateChanged(int)));

		connect(le_scale, SIGNAL(editingFinished()), SLOT(le_scale_editingFinished()));
		connect(le_scale, SIGNAL(returnPressed()), SLOT(le_scale_editingFinished()));

		connect(sl_thickness, SIGNAL(valueChanged(int)), SLOT(sl_thickness_valueChanged(int)));
		connect(sl_detector, SIGNAL(valueChanged(int)), SLOT(sl_detector_valueChanged(int)));
		connect(ckb_show_coh, SIGNAL(stateChanged(int)), SLOT(ckb_show_coh_stateChanged(int)));

		connect(le_rep_x, SIGNAL(textChanged(QString)), SLOT(le_rep_xy_textChanged(QString)));
		connect(le_rep_y, SIGNAL(textChanged(QString)), SLOT(le_rep_xy_textChanged(QString)));

		connect(image, SIGNAL(sg_save_as_binary()), this, SLOT(save_as_binary()));
		/******************************************************/
		draw_current_settings();
	}

public slots:
	void save_as_binary()
	{
		QString FileName = QFileDialog::getSaveFileName(0,"Save file", QDir::currentPath(),
		"Binary file (*.bin );Binary (*.bin)",
		new QString("Binary file (*.bin )"));

		if (!FileName.isEmpty())
		{
		   //mt::Grid_2d<float> grid(output_multislice.nx, output_multislice.grid_2d.ny);
		   //mt::write_mat_binary_matrix(FileName, grid, data);
		}
	}

	void set_normalSize()
	{
		adjustSize();
	}

	void qt_image_mouseMoveEvent(QPoint Point)
	{
//		lb_Pos->setText(QString(" Position: X = %1, Y = %2").arg(Point.rx()).arg(Point.ry()));
//		lb_Rel_Pos->setText(QString(" Position: X = %1, Y = %2").arg(Point.rx()-MouSetPos.rx()).arg(Point.ry()-MouSetPos.ry()));
	}

	void qt_image_mousePressEvent(QPoint Point)
	{
//		lb_Rel_Ori->setText(QString(" Position: X = %1, Y = %2").arg(Point.rx()).arg(Point.ry()));
//		MouSetPos = Point;
	}

	void qt_image_mouseReleaseEvent(QPoint Point)
	{

	}

	void cb_colormap_currentIndexChanged(int index)
	{
		auto show_cdata = QVar_to_enum<mt::eShow_CData>(cb_show_ctype->currentData());
		auto colormap = QVar_to_enum<eColormap>(cb_colormap->itemData(index));
		auto scale_type = QVar_to_enum<eScale_Type>(cb_scale->currentData());
		double scale_factor = le_scale->text().toDouble();

		int ithk = sl_thickness->value()-1;
		int idet = sl_detector->value()-1;

		draw(show_cdata, colormap, scale_type, scale_factor, ithk, idet);
	}

	void cb_show_ctype_currentIndexChanged(int index)
	{
		auto show_cdata = QVar_to_enum<mt::eShow_CData>(cb_show_ctype->itemData(index));
		auto colormap = QVar_to_enum<eColormap>(cb_colormap->currentData());

		if(!(show_cdata==mt::eSCD_CMod))
		{
			cb_scale->setCurrentIndex(0);
		}
		auto scale_type = QVar_to_enum<eScale_Type>(cb_scale->currentData());
		double scale_factor = le_scale->text().toDouble();

		int ithk = sl_thickness->value()-1;
		int idet = sl_detector->value()-1;

		draw(show_cdata, colormap, scale_type, scale_factor, ithk, idet);
	}

	void cb_scale_currentIndexChanged(int index)
	{
		auto show_cdata = QVar_to_enum<mt::eShow_CData>(cb_show_ctype->currentData());
		auto colormap = QVar_to_enum<eColormap>(cb_colormap->currentData());
		auto scale_type = QVar_to_enum<eScale_Type>(cb_scale->itemData(index));
		double scale_factor = le_scale->text().toDouble();

		auto b_st = (scale_type==eSLT_Power);
		le_scale->setVisible(b_st);

		int ithk = sl_thickness->value()-1;
		int idet = sl_detector->value()-1;

		draw(show_cdata, colormap, scale_type, scale_factor, ithk, idet);
	}

	void le_scale_editingFinished()
	{
		auto show_cdata = QVar_to_enum<mt::eShow_CData>(cb_show_ctype->currentData());
		auto colormap = QVar_to_enum<eColormap>(cb_colormap->currentData());
		auto scale_type = QVar_to_enum<eScale_Type>(cb_scale->currentData());
		double scale_factor = le_scale->text().toDouble();

		int ithk = sl_thickness->value()-1;
		int idet = sl_detector->value()-1;

		draw(show_cdata, colormap, scale_type, scale_factor, ithk, idet);
	}

	void le_rep_xy_textChanged(const QString &)
	{
		auto show_cdata = QVar_to_enum<mt::eShow_CData>(cb_show_ctype->currentData());
		auto colormap = QVar_to_enum<eColormap>(cb_colormap->currentData());
		auto scale_type = QVar_to_enum<eScale_Type>(cb_scale->currentData());
		double scale_factor = le_scale->text().toDouble();

		int ithk = sl_thickness->value()-1;
		int idet = sl_detector->value()-1;

		int rep_nx = le_rep_x->text().toInt();
		int rep_ny = le_rep_y->text().toInt();

		auto pmo = (ckb_show_coh->isChecked())?mt::eFMO_Coherent:mt::eFMO_Total;
		data = output_multislice.extract_data(pmo, show_cdata, ithk, idet);
		data = replicate(rep_nx, rep_ny, data);

		image->set_data_size(rep_nx*output_multislice.nx, rep_ny*output_multislice.ny);
		image->draw(&data, colormap, scale_type, scale_factor);
	}

	void ckb_block_g0_stateChanged(int state)
	{
		draw_current_settings();
	}

	void ckb_show_coh_stateChanged(int state)
	{
		if(state==Qt::Checked)
		{
			if(output_multislice.is_ot_m2psi_tot_psi_coh())
			{
				lb_show_ctype->setVisible(true);
				cb_show_ctype->setVisible(true);
			}
		}
		else
		{
			auto bb_complex = is_output_complex_alone();
			lb_show_ctype->setVisible(bb_complex);
			cb_show_ctype->setVisible(bb_complex);
		}
		draw_current_settings();
	}

	void sl_thickness_valueChanged(int value)
	{
		if(value>output_multislice.thick.size()) return;

		auto show_cdata = QVar_to_enum<mt::eShow_CData>(cb_show_ctype->currentData());
		auto colormap = QVar_to_enum<eColormap>(cb_colormap->currentData());
		auto scale_type = QVar_to_enum<eScale_Type>(cb_scale->currentData());
		double scale_factor = le_scale->text().toDouble();

		int ithk = value-1;
		int idet = sl_detector->value()-1;

		lb_thickness_value->setText(thk_value_qstr(output_multislice.thick[ithk]));

		draw(show_cdata, colormap, scale_type, scale_factor, ithk, idet);
	}

	void sl_detector_valueChanged(int value)
	{
		if(value>output_multislice.ndetector) return;

		auto show_cdata = QVar_to_enum<mt::eShow_CData>(cb_show_ctype->currentData());
		auto colormap = QVar_to_enum<eColormap>(cb_colormap->currentData());
		auto scale_type = QVar_to_enum<eScale_Type>(cb_scale->currentData());
		double scale_factor = le_scale->text().toDouble();

		int ithk = sl_thickness->value()-1;
		int idet = value-1;

		auto lambda = mt::get_lambda(output_multislice.E_0);
		auto g_inner = output_multislice.detector.g_inner[idet];
		auto g_outer = output_multislice.detector.g_outer[idet];

		lb_detector_value->setText(det_value_qstr(lambda, g_inner, g_outer));

		int rep_nx = le_rep_x->text().toInt();
		int rep_ny = le_rep_y->text().toInt();

		auto pmo = (ckb_show_coh->isChecked())?mt::eFMO_Coherent:mt::eFMO_Total;
		data = output_multislice.extract_data(pmo, show_cdata, ithk, idet);
		data = replicate(rep_nx, rep_ny, data);

		image->set_data_size(rep_nx*output_multislice.nx, rep_ny*output_multislice.ny);
		image->draw(&data, colormap, scale_type, scale_factor);
	}

	void draw_current_settings()
	{
		auto show_cdata = QVar_to_enum<mt::eShow_CData>(cb_show_ctype->currentData());
		auto colormap = QVar_to_enum<eColormap>(cb_colormap->currentData());
		auto scale_type = QVar_to_enum<eScale_Type>(cb_scale->currentData());
		double scale_factor = le_scale->text().toDouble();

		int ithk = sl_thickness->value()-1;
		int idet = sl_detector->value()-1;

		draw(show_cdata, colormap, scale_type, scale_factor, ithk, idet);
	}

private:
	/****************************************************/
	QLabel *lb_colormap;
	QComboBox *cb_colormap;

	QCheckBox *ckb_show_coh;

	QLabel *lb_scale;
	QComboBox *cb_scale;
	QLineEdit *le_scale;

	QLabel *lb_show_ctype;
	QComboBox *cb_show_ctype;

	QCheckBox *ckb_block_g0;

	QLabel *lb_thickness;
	QSlider *sl_thickness;
	QLabel *lb_thickness_value;

	QLabel *lb_detector;
	QSlider *sl_detector;
	QLabel *lb_detector_value;

	/****************************************************/
	QLabel *lb_rep_x;
	QLineEdit *le_rep_x;
	QLabel *lb_rep_y;
	QLineEdit *le_rep_y;

	QLabel *lb_source_size;
	QLineEdit *le_source_size;

	QLabel *lb_pos_pixel;
	QLabel *lb_pos_real;

	/****************************************************/
	Data_Image *image;

	QPoint MouSetPos;

	QSize window_size;

	/****************************************************/
	mt::Output_Multislice<float> output_multislice;

	host_vector<float> data;

	void set_title()
	{
		auto sim_type = output_multislice.simulation_type;

		switch(sim_type)
		{
			case mt::eTEMST_STEM:
				setWindowTitle(tr("Scanning transmission electron microscopy"));
			break;
			case mt::eTEMST_ISTEM:
				setWindowTitle(tr("Imaging scanning transmission electron microscopy"));
			break;
			case mt::eTEMST_CBED:
				setWindowTitle(tr("Convergent beam electron diffraction"));
			break;
			case mt::eTEMST_CBEI:
				setWindowTitle(tr("Convergent beam electron imaging"));
			break;
			case mt::eTEMST_HRTEM:
				setWindowTitle(tr("High resolution transmission electron microscopy"));
			break;
			case mt::eTEMST_ED:
				setWindowTitle(tr("Electron diffraction"));
			break;
			case mt::eTEMST_PED:
				setWindowTitle(tr("Precession electron diffraction"));
			break;
			case mt::eTEMST_HCTEM:
				setWindowTitle(tr("Hollow cone TEM"));
			break;
			case mt::eTEMST_EWFS:
				setWindowTitle(tr("Exit wave in Fourier space"));
			break;
			case mt::eTEMST_EWRS:
				setWindowTitle(tr("Exit wave in real space"));
			break;
			case mt::eTEMST_EELS:
				setWindowTitle(tr("STEM electron energy loss spectroscopy"));
			break;
			case mt::eTEMST_EFTEM:
				setWindowTitle(tr("Energy filtered TEM"));
			break;
			case mt::eTEMST_IWFS:
				setWindowTitle(tr("Incident wave in real space"));
			break;
			case mt::eTEMST_IWRS:
				setWindowTitle(tr("Incident wave in Fourier space"));
			break;
			case mt::eTEMST_PPRS:
				setWindowTitle(tr("Projected potential in real space"));
			break;
			case mt::eTEMST_TFRS:
				setWindowTitle(tr("Transmission function in real space"));
			break;
		}
	}

	void create_wigets()
	{
		lb_colormap = new QLabel(tr("Colormap"));

		cb_colormap = new QComboBox();
		cb_colormap->addItem("Gray", QVariant(eCM_Gray));
		cb_colormap->addItem("Hot", QVariant(eCM_Hot));
		cb_colormap->addItem("Cool", QVariant(eCM_Cool));
		cb_colormap->addItem("Jet", QVariant(eCM_Jet));
		cb_colormap->addItem("Copper", QVariant(eCM_Copper));
		cb_colormap->addItem("Summer", QVariant(eCM_Summer));
		cb_colormap->addItem("Autumn", QVariant(eCM_Autumn));
		cb_colormap->addItem("Winter", QVariant(eCM_Winter));

		/***************************************************/
		ckb_show_coh = new QCheckBox(tr("Show Coherent wave"));
		ckb_show_coh->setStatusTip(tr("Show Coherent wave"));

		/***************************************************/
		lb_scale = new QLabel(tr("Scale"));

		cb_scale = new QComboBox();
		cb_scale->addItem("Linear", QVariant(eSLT_Linear));
		cb_scale->addItem("Power", QVariant(eSLT_Power));
		cb_scale->addItem("Log", QVariant(eSLT_Log));

		le_scale = new QLineEdit;
		le_scale->setStatusTip(tr("Scaling factor [0-1]"));
		le_scale->setMaximumWidth(70);
		le_scale->setValidator(new QDoubleValidator(0, 1.0, 3));
		le_scale->setText("0.500");

		/***************************************************/
		lb_show_ctype = new QLabel(tr("Complex"));

		cb_show_ctype = new QComboBox();
		cb_show_ctype->addItem("Module", QVariant(mt::eSCD_CMod));
		cb_show_ctype->addItem("Phase", QVariant(mt::eSCD_CPhase));
		cb_show_ctype->addItem("Real", QVariant(mt::eSCD_CReal));
		cb_show_ctype->addItem("Imag", QVariant(mt::eSCD_CImag));

		/***************************************************/
		ckb_block_g0 = new QCheckBox(tr("Block g(0))"));
		ckb_block_g0->setStatusTip(tr("Block g(0)"));

		/***************************************************/
		lb_thickness = new QLabel(tr(" Thickness "));
		lb_thickness->setFrameStyle(QFrame::Box | QFrame::Sunken);
		lb_thickness->setAlignment(Qt::AlignCenter);

		sl_thickness = new QSlider(Qt::Horizontal);
		sl_thickness->setMinimumWidth(25);

		lb_thickness_value = new QLabel(thk_value_qstr(0));
		lb_thickness_value->setFrameStyle(QFrame::Box | QFrame::Sunken);
		lb_thickness_value->setAlignment(Qt::AlignCenter);

		/***************************************************/
		lb_detector = new QLabel(tr(" Detector "));
		lb_detector->setFrameStyle(QFrame::Box | QFrame::Sunken);
		lb_detector->setAlignment(Qt::AlignCenter);

		sl_detector = new QSlider(Qt::Horizontal);
		sl_detector->setMinimumWidth(25);

		lb_detector_value = new QLabel(det_value_qstr(0.019, 0, 0));
		lb_detector_value->setFrameStyle(QFrame::Box | QFrame::Sunken);
		lb_detector_value->setAlignment(Qt::AlignCenter);

		/***************************************************/
		/***************************************************/
		auto lyh_top = new QHBoxLayout;
		lyh_top->setContentsMargins(0, 0, 0, 0);
		lyh_top->setSizeConstraint(QLayout::SetMinimumSize);
		lyh_top->addWidget(lb_colormap);
		lyh_top->addWidget(cb_colormap);
		lyh_top->addWidget(lb_scale);
		lyh_top->addWidget(cb_scale);
		lyh_top->addWidget(le_scale);
		lyh_top->addWidget(ckb_show_coh);
		lyh_top->addWidget(lb_show_ctype);
		lyh_top->addWidget(cb_show_ctype);
		lyh_top->addWidget(ckb_block_g0);
		lyh_top->addWidget(lb_thickness);
		lyh_top->addWidget(sl_thickness);
		lyh_top->addWidget(lb_thickness_value);
		lyh_top->addWidget(lb_detector);
		lyh_top->addWidget(sl_detector);
		lyh_top->addWidget(lb_detector_value);
		lyh_top->addStretch();

		auto gb_top = new QGroupBox;
		gb_top->setMouseTracking(true);
		gb_top->setLayout(lyh_top);
		gb_top->setStyleSheet("QGroupBox{border: 1px solid gray; \n border-radius: 2px; \n margin-top: 0px;}");
		gb_top->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);

		/***************************************************/
		lb_rep_x = new QLabel("Rep. x");

		le_rep_x = new QLineEdit;
		le_rep_x->setStatusTip(tr("Replicate along x"));
		le_rep_x->setMaximumWidth(70);
		le_rep_x->setValidator(new QIntValidator(0, 50));
		le_rep_x->setText("1");

		/***************************************************/
		lb_rep_y = new QLabel("Rep. y");

		le_rep_y = new QLineEdit;
		le_rep_y->setStatusTip(tr("Replicate along x"));
		le_rep_y->setMaximumWidth(70);
		le_rep_y->setValidator(new QIntValidator(0, 50));
		le_rep_y->setText("1");

		/***************************************************/
		lb_source_size = new QLabel("Source Size [Å]");

		le_source_size = new QLineEdit;
		le_source_size->setStatusTip(tr("Source size (FWHM)"));
		le_source_size->setMaximumWidth(70);
		le_source_size->setValidator(new QDoubleValidator(0, 100, 3));
		le_source_size->setText("0.800");
		le_source_size->setVisible(false);

		/***************************************************/
		lb_pos_pixel = new QLabel(" (0, 0) Pixels");
		lb_pos_real = new QLabel(" (0, 0) [Å]");

		/***************************************************/
		/***************************************************/
		auto lyh_bottom = new QHBoxLayout;
		lyh_bottom->setContentsMargins(0, 0, 0, 0);
		lyh_bottom->setSizeConstraint(QLayout::SetMinimumSize);
		lyh_bottom->addWidget(lb_rep_x);
		lyh_bottom->addWidget(le_rep_x);
		lyh_bottom->addWidget(lb_rep_y);
		lyh_bottom->addWidget(le_rep_y);
		lyh_bottom->addWidget(lb_source_size);
		lyh_bottom->addWidget(le_source_size);
		lyh_bottom->addStretch();
		lyh_bottom->addWidget(lb_pos_pixel);
		lyh_bottom->addWidget(lb_pos_real);

		auto gb_bottom = new QGroupBox;
		gb_bottom->setMouseTracking(true);
		gb_bottom->setLayout(lyh_bottom);
		gb_bottom->setStyleSheet("QGroupBox{border: 1px solid gray; \n border-radius: 2px; \n margin-top: 0px;}");
		gb_bottom->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);

		/***************************************************/
		image = new Data_Image;
		image->setMouseTracking(true);
		image->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

		/***************************************************/
		auto lyv_main = new QVBoxLayout;
		lyv_main->setSpacing(0);
		lyv_main->setMargin(0);
		lyv_main->setContentsMargins(0, 0, 0, 0);
		lyv_main->addWidget(gb_top);
		lyv_main->addWidget(image);
		lyv_main->addWidget(gb_bottom);
		setLayout(lyv_main);

//		connect(q_image, SIGNAL(Signal_NormalSize(QSize)), this, SLOT(qt_image_NormalSize(QSize)));
//		connect(q_image, SIGNAL(Signal_mouseMoveEvent(QPoint)), this, SLOT(qt_image_mouseMoveEvent(QPoint)));
//		connect(q_image, SIGNAL(Signal_mousePressEvent(QPoint)), this, SLOT(qt_image_mousePressEvent(QPoint)));
//		connect(q_image, SIGNAL(Signal_mouseReleaseEvent(QPoint)), this, SLOT(qt_image_mouseReleaseEvent(QPoint)));
//		connect(cb_colormap, SIGNAL(currentIndexChanged(int)), this, SLOT(on_cb_colormap_currentIndexChanged(int)));
	}

	void activate_options()
	{
		auto b_complex = is_output_complex_alone();

		lb_show_ctype->setVisible(b_complex);
		cb_show_ctype->setVisible(b_complex);

		auto b_pwr = output_multislice.is_CBED_ED_EWFS_PED();
		le_scale->setVisible(b_pwr);

		if(b_pwr)
		{
			cb_scale->setCurrentIndex(1);
			le_scale->setText(QString::number(0.25, 'f', 3));
		}
		else
		{
			cb_scale->setCurrentIndex(0);
		}

		ckb_show_coh->setVisible(is_two_ouputs_present());

		auto b_g0 = (output_multislice.is_ED() || output_multislice.is_EWFS() || output_multislice.is_PED()) && output_multislice.is_plane_wave();

		ckb_block_g0->setChecked(true);
		ckb_block_g0->setVisible(b_g0);

		auto b_thk_sl = (output_multislice.thick.size()>1);
		sl_thickness->setVisible(b_thk_sl);

		auto b_thk_val = is_thickness_available();
		lb_thickness->setVisible(b_thk_val);
		lb_thickness_value->setVisible(b_thk_val);

		/**************************************************************/
		auto b_stem = output_multislice.is_STEM();

		lb_rep_x->setVisible(b_stem);
		le_rep_x->setVisible(b_stem);
		lb_rep_y->setVisible(b_stem);
		le_rep_y->setVisible(b_stem);
		lb_source_size->setVisible(false);
		le_source_size->setVisible(false);

		auto b_det_sl = b_stem && (output_multislice.ndetector>1);
		sl_detector->setVisible(b_det_sl);

		auto b_det_val = b_stem;
		lb_detector->setVisible(b_det_val);
		lb_detector_value->setVisible(b_det_val);
	}

	bool is_two_ouputs_present()
	{
		bool b = output_multislice.is_ot_image_tot_coh()||output_multislice.is_ot_m2psi_tot_coh();
		b = b || output_multislice.is_ot_m2psi_tot_psi_coh();

		return b;
	}

	bool is_output_complex_alone()
	{
		bool b = output_multislice.is_ot_psi_coh()||output_multislice.is_ot_psi_0();
		b = b || output_multislice.is_ot_trans();

		return b;
	}

	bool is_thickness_available()
	{
		bool bb =  output_multislice.is_IWFS_IWRS() || output_multislice.is_PPFS_PPRS();
		bb = bb || output_multislice.is_TFFS_TFRS() || output_multislice.is_PropFS_PropRS();

		return !bb;
	}

	host_vector<float> replicate(int nx_r, int ny_r, host_vector<float> &data)
	{
		int nx = output_multislice.nx;
		int ny = output_multislice.ny;

		int nx_t = nx_r*nx;
		int ny_t = ny_r*ny;

		host_vector<float> data_t(nx_t*ny_t);

		int ix_t = 0;
		for(auto ix_r=0; ix_r<nx_r; ix_r++)
		{
			for(auto ix=0; ix<nx; ix++)
			{
				int iy_t = 0;
				for(auto iy_r=0; iy_r<ny_r; iy_r++)
				{
					for(auto iy=0; iy<ny; iy++)
					{
						data_t[ix_t*ny_t+iy_t] = data[ix*ny+iy];
						iy_t++;
					}
				}
				ix_t++;
			}
		}
		return data_t;
	}

	void draw(mt::eShow_CData show_cdata, eColormap colormap,
	eScale_Type scale_type, double scale_factor, int ithk, int idet=0)
	{
		auto pmo = (ckb_show_coh->isChecked())?mt::eFMO_Coherent:mt::eFMO_Total;
		data = output_multislice.extract_data(pmo, show_cdata, ithk, idet);


		auto b_g0 = (output_multislice.is_ED() || output_multislice.is_EWFS() || output_multislice.is_PED()) && output_multislice.is_plane_wave();
		if(b_g0 && ckb_block_g0->isChecked())
		{
			int ix = output_multislice.nx/2;
			int iy = output_multislice.ny/2;
			int ny = output_multislice.ny;

			data[ix*ny+iy] = 0;
		}

		image->draw(&data, colormap, scale_type, scale_factor);
	}

//protected:
//	void closeEvent(QCloseEvent *event);

//signals:
//	void Signal_CloseWindow(QString name);

};

#endif // DATA_VIEWER_H
