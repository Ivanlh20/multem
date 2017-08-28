#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "types.cuh"
#include "atomic_data.hpp"
#include "atom_data.hpp"
#include "input_multislice.cuh"
#include "output_multislice.hpp"

#include <QMainWindow>
#include <QThread>
#include <QTimer>

QT_BEGIN_NAMESPACE
class QAction;
class QDialogButtonBox;
class QGroupBox;
class QFrame;
class QLabel;
class QLineEdit;
class QMenu;
class QMenuBar;
class QPushButton;
class QTextEdit;
class QRadioButton;
class QComboBox;
class QSpinBox;
class QCheckBox;
class QTabWidget;
class QProgressBar;
class MGroupBox_V;
class MGroupBox_H;
class MGroupBox_G;
class MLine_H;
class MGroupBox_HE;
class MGroupBox_VE;
class MGroupBox_GE;
class MFrame_V;
class MFrame_H;
class MFrame_G;

QT_END_NAMESPACE

class RS_Thread;
class PB_Timer;

/************************************************/
class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow();

	RS_Thread *rs_thread;
	PB_Timer *pb_timer;

	void run_sim();

//public signals:
//	void close_childs();

private slots:
	void open_sim();
	void save_sim();
	void save_as_sim();
	void exit();
	void about();
	void default_sim();
	void start_sim();
	void stop_sim();

public slots:
	void progress_run_update();
	void show_data();

	void cb_device_currentIndexChanged(int index);
	void cb_elec_spec_int_model_currentIndexChanged(int index);

	void pb_spec_load_released();

	void le_spec_lx_ly_textChanged(const QString &);
	void pb_spec_recenter_released();
	void cb_spec_rot_center_type_currentIndexChanged(int index);
	void cb_pot_slic_type_currentIndexChanged(int index);
	void cb_thk_type_currentIndexChanged(int index);
	void le_pot_slic_thick_editingFinished();

	void cb_sim_type_currentIndexChanged(int index);

	void cb_iw_type_currentIndexChanged(int index);
	void pb_iw_load_released();

	void cb_stem_det_type_currentIndexChanged(int index);
	void le_stem_det_n_det_editingFinished();
	void cb_stem_det_k_det_currentIndexChanged(int index);
	void le_stem_det_ideal_editingFinished();
	void pb_stem_det_exp_load_released();

	void cb_eels_element_currentIndexChanged(int index);

	void cb_illu_model_currentIndexChanged(int index);
	void cb_cl_sa_c_10_zero_currentIndexChanged(int index);
	void cb_ol_sa_c_10_zero_currentIndexChanged(int index);
	void pb_cl_sa_c_10_opt_released();
	void pb_ol_sa_c_10_opt_released();

private:
	void create_menubar();

	QLabel *central_widget;
	QMenu *file_menu;
	QAction *act_open_sim;
	QAction *act_save_sim;
	QAction *act_save_as_sim;
	QAction *act_exit;
	QAction *act_default_sim;
	QAction *act_start_sim;
	QAction *act_stop_sim;

	QProgressBar *pb_progress_run;
	/****************General*****************/
	MGroupBox_G *gbg_device;
	QComboBox *cb_device;

	QLabel *lb_precision;
	QComboBox *cb_precision;

	QLabel *lb_gpu_card;
	QComboBox *cb_gpu_card;

	QLabel *lb_nthreads;
	QSpinBox *sb_nthreads;

	MGroupBox_G *gbg_elec_spec_int_model;
	QComboBox *cb_elec_spec_int_model;
	QLabel *lb_potential_type;
	QComboBox *cb_potential_type;

	MGroupBox_G *gbg_elec_phonon_int_model;
	QComboBox *cb_pn_model;
	QCheckBox *ckb_pn_coh_contrib;
	QLabel *lb_pn_dim;
	QComboBox *cb_pn_dim;
	QLabel *lb_pn_nconf;
	QLineEdit *le_pn_nconf;
	QCheckBox *ckb_pn_single_conf;
	QLabel *lb_pn_seed;
	QLineEdit *le_pn_seed;
	QCheckBox *ckb_pn_auto_seed;
	QWidget *wg_electron_phonon;

	QDockWidget *dw_general;

	/***********************Specimen**********************/
	QLabel *lb_spec_file;
	QLineEdit *le_spec_file;
	QPushButton *pb_spec_load;
	QPushButton *pb_spec_show;

	MGroupBox_G *gbg_spec_info;

	QLabel *lb_spec_n_atom_types;
	QLabel *lb_spec_lx;
	QLineEdit *le_spec_lx;
	QLabel *lb_spec_ly;
	QLineEdit *le_spec_ly;
	QPushButton *pb_spec_recenter;

	MGroupBox_G *gbg_spec_rot;

	QLabel *lb_spec_rot_theta;
	QLineEdit *le_spec_rot_theta;
	QLabel *lb_spec_rot_u;
	QLineEdit *le_spec_rot_ux;
	QLineEdit *le_spec_rot_uy;
	QLineEdit *le_spec_rot_uz;
	QLabel *lb_spec_rot_center_type;
	QComboBox *cb_spec_rot_center_type;
	QLabel *lb_spec_rot_center_p;
	QLineEdit *le_spec_rot_center_px;
	QLineEdit *le_spec_rot_center_py;
	QLineEdit *le_spec_rot_center_pz;

	/*****************************************************/
	MGroupBox_G *gbg_thk;
	QLabel *lb_thk_type;
	QComboBox *cb_thk_type;

	QLabel *lb_thk_0;
	QLineEdit *le_thk_0;
	QLabel *lb_thk_d;
	QLineEdit *le_thk_d;
	QLabel *lb_thk_e;
	QLineEdit *le_thk_e;

	MGroupBox_G *gbg_potential_slic;

	QComboBox *cb_pot_slic_type;
	QLabel *lb_pot_slic_thick;
	QLineEdit *le_pot_slic_thick;
	QPushButton *pb_pot_slic_view;

	MGroupBox_G *gbg_samp;

	QLabel *lb_samp_nx;
	QLineEdit *le_samp_nx;
	QLabel *lb_samp_ny;
	QLineEdit *le_samp_ny;
	QCheckBox *ckb_samp_bandwidth;

	QDockWidget *dw_specimen;

	/****************Incident wave*****************/
	QLabel *lb_iw_type;
	QComboBox *cb_iw_type;

	QLabel *lb_iw_file;
	QLineEdit *le_iw_file;
	QPushButton *pb_iw_load;

	QLabel *lb_iw_x;
	QLineEdit *le_iw_x;
	QLabel *lb_iw_y;
	QLineEdit *le_iw_y;

	QPushButton *pb_iw_p;

	MGroupBox_G *gbg_iw;

	QWidget *wg_iw;

	/*********************stem**********************/
	QLabel *lb_stem_sc_type;
	QComboBox *cb_stem_sc_type;
	QCheckBox *ckb_stem_sc_incl_last_point;

	QLabel *lb_stem_sc_n_points;
	QLineEdit *le_stem_sc_n_points;
	QLabel *lb_stem_sc_np_assign;
	QComboBox *cb_stem_sc_np_assign;

	QLabel *lb_stem_sc_px_0;
	QLineEdit *le_stem_sc_px_0;
	QLabel *lb_stem_sc_py_0;
	QLineEdit *le_stem_sc_py_0;
	QPushButton *pb_stem_sc_p_0;

	QLabel *lb_stem_sc_px_e;
	QLineEdit *le_stem_sc_px_e;
	QLabel *lb_stem_sc_py_e;
	QLineEdit *le_stem_sc_py_e;
	QPushButton *pb_stem_sc_p_e;

	MGroupBox_G *gbg_stem_scanning;

	QComboBox *cb_stem_det_type;
	QLabel *lb_stem_det_n_det;
	QLineEdit *le_stem_det_n_det;
	QLabel *lb_stem_det_k_det;
	QComboBox *cb_stem_det_k_det;
	QPushButton *pb_stem_det_view;

	MGroupBox_G *gbg_stem_detector;

	QFrame *fr_stem_detector_ideal;
	QLabel *lb_stem_det_ideal_inner_angle;
	QLineEdit *le_stem_det_ideal_inner_angle;
	QLabel *lb_stem_det_ideal_outer_angle;
	QLineEdit *le_stem_det_ideal_outer_angle;

	QFrame *fr_stem_detector_exp;
	QLabel *lb_stem_det_exp_file;
	QLineEdit *le_stem_det_exp_file;
	QPushButton *pb_stem_det_exp_load;

	QWidget *wg_stem;

	/*************PED/HCTEM**************/
	QLabel *lb_pcs_tilt_angle;
	QLineEdit *le_pcs_tilt_angle;

	QLabel *lb_pcs_n_rotation;
	QLineEdit *le_pcs_n_rotation;

	MGroupBox_G *gbg_pcs;

	QWidget *wg_pcs;

	/*********************EFTEM/EELS***********************/
	QLabel *lb_eels_element;
	QComboBox *cb_eels_element;

	QLabel *lb_eels_energy;
	QComboBox *cb_eels_energy;

	QLabel *lb_eels_coll_angle;
	QLineEdit *le_eels_coll_angle;

	QLabel *lb_eels_m_selection;
	QComboBox *cb_eels_m_selection;

	QLabel *lb_eels_channelling_type;
	QComboBox *cb_eels_channelling_type;

	MGroupBox_G *gbg_eels;

	QWidget *wg_eels;

	/*************************PP/TF************************/
	QLabel *lb_pptf_n_slices;

	QLabel *lb_pptf_slice;
	QComboBox *cb_pptf_slice;

	MGroupBox_H *gbh_pptf;

	QWidget *wg_pptf;

	/*******************************************************/
	MGroupBox_G *gbg_simulation;

	QLabel *lb_sim_type;
	QComboBox *cb_sim_type;
	QDockWidget *dw_simulation;

	/************Microscope parameters**************/
	QLabel *lb_mic_name;
	QComboBox *cb_mic_name;
	QLabel *lb_mic_acc_vol;
	QLineEdit *le_mic_acc_vol;
	QPushButton *pb_mic_view;

	/**********************************************/
	QLabel *lb_tilt_theta;
	QLineEdit *le_tilt_theta;
	QLabel *lb_tilt_phi;
	QLineEdit *le_tilt_phi;

	/**********************************************/
	QLabel *lb_illu_model;
	QComboBox *cb_illu_model;

	QLabel *lb_illu_beta;
	QLineEdit *le_illu_beta;
	QLabel *lb_illu_nbeta;
	QLineEdit *le_illu_nbeta;

	MGroupBox_G *gbg_illu;

	/*************Condenser lens***************/
	QLabel *lb_cl_sa_c_30;
	QLineEdit *le_cl_sa_c_30;
	QLabel *lb_cl_sa_c_50;
	QLineEdit *le_cl_sa_c_50;
	QLabel *lb_cl_sa_m;
	QLineEdit *le_cl_sa_m;

	QLabel *lb_cl_sa_c_10;
	QLineEdit *le_cl_sa_c_10;
	QPushButton *pb_cl_sa_c_10_opt;
	QLabel *lb_cl_sa_c_10_zero;
	QComboBox *cb_cl_sa_c_10_zero;
	QLabel *lb_cl_sa_c_10_z;
	QLineEdit *le_cl_sa_c_10_z;

	QLabel *lb_cl_sa_delta;
	QLineEdit *le_cl_sa_delta;
	QLabel *lb_cl_sa_ndelta;
	QComboBox *cb_cl_sa_ndelta;

	MGroupBox_GE *gbge_cl_nsa;

	QLabel *lb_cl_c_12;
	QLineEdit *le_cl_c_12;
	QLabel *lb_cl_phi_12;
	QLineEdit *le_cl_phi_12;

	QLabel *lb_cl_c_21;
	QLineEdit *le_cl_c_21;
	QLabel *lb_cl_phi_21;
	QLineEdit *le_cl_phi_21;
	QLabel *lb_cl_c_23;
	QLineEdit *le_cl_c_23;
	QLabel *lb_cl_phi_23;
	QLineEdit *le_cl_phi_23;

	QLabel *lb_cl_c_32;
	QLineEdit *le_cl_c_32;
	QLabel *lb_cl_phi_32;
	QLineEdit *le_cl_phi_32;
	QLabel *lb_cl_c_34;
	QLineEdit *le_cl_c_34;
	QLabel *lb_cl_phi_34;
	QLineEdit *le_cl_phi_34;

	QLabel *lb_cl_c_41;
	QLineEdit *le_cl_c_41;
	QLabel *lb_cl_phi_41;
	QLineEdit *le_cl_phi_41;
	QLabel *lb_cl_c_43;
	QLineEdit *le_cl_c_43;
	QLabel *lb_cl_phi_43;
	QLineEdit *le_cl_phi_43;
	QLabel *lb_cl_c_45;
	QLineEdit *le_cl_c_45;
	QLabel *lb_cl_phi_45;
	QLineEdit *le_cl_phi_45;

	QLabel *lb_cl_c_52;
	QLineEdit *le_cl_c_52;
	QLabel *lb_cl_phi_52;
	QLineEdit *le_cl_phi_52;
	QLabel *lb_cl_c_54;
	QLineEdit *le_cl_c_54;
	QLabel *lb_cl_phi_54;
	QLineEdit *le_cl_phi_54;
	QLabel *lb_cl_c_56;
	QLineEdit *le_cl_c_56;
	QLabel *lb_cl_phi_56;
	QLineEdit *le_cl_phi_56;

	MGroupBox_H *gbh_cl_ar;
	QLabel *lb_cl_ar_min;
	QLineEdit *le_cl_ar_min;
	QLabel *lb_cl_ar_max;
	QLineEdit *le_cl_ar_max;

	/*************Objective lens***************/
	QLabel *lb_ol_sa_c_30;
	QLineEdit *le_ol_sa_c_30;
	QLabel *lb_ol_sa_c_50;
	QLineEdit *le_ol_sa_c_50;
	QLabel *lb_ol_sa_m;
	QLineEdit *le_ol_sa_m;

	QLabel *lb_ol_sa_c_10;
	QLineEdit *le_ol_sa_c_10;
	QPushButton *pb_ol_sa_c_10_opt;
	QLabel *lb_ol_sa_c_10_zero;
	QComboBox *cb_ol_sa_c_10_zero;
	QLabel *lb_ol_sa_c_10_z;
	QLineEdit *le_ol_sa_c_10_z;

	QLabel *lb_ol_sa_delta;
	QLineEdit *le_ol_sa_delta;
	QLabel *lb_ol_sa_ndelta;
	QComboBox *cb_ol_sa_ndelta;

	MGroupBox_GE *gbge_ol_nsa;

	QLabel *lb_ol_c_12;
	QLineEdit *le_ol_c_12;
	QLabel *lb_ol_phi_12;
	QLineEdit *le_ol_phi_12;

	QLabel *lb_ol_c_21;
	QLineEdit *le_ol_c_21;
	QLabel *lb_ol_phi_21;
	QLineEdit *le_ol_phi_21;
	QLabel *lb_ol_c_23;
	QLineEdit *le_ol_c_23;
	QLabel *lb_ol_phi_23;
	QLineEdit *le_ol_phi_23;

	QLabel *lb_ol_c_32;
	QLineEdit *le_ol_c_32;
	QLabel *lb_ol_phi_32;
	QLineEdit *le_ol_phi_32;
	QLabel *lb_ol_c_34;
	QLineEdit *le_ol_c_34;
	QLabel *lb_ol_phi_34;
	QLineEdit *le_ol_phi_34;

	QLabel *lb_ol_c_41;
	QLineEdit *le_ol_c_41;
	QLabel *lb_ol_phi_41;
	QLineEdit *le_ol_phi_41;
	QLabel *lb_ol_c_43;
	QLineEdit *le_ol_c_43;
	QLabel *lb_ol_phi_43;
	QLineEdit *le_ol_phi_43;
	QLabel *lb_ol_c_45;
	QLineEdit *le_ol_c_45;
	QLabel *lb_ol_phi_45;
	QLineEdit *le_ol_phi_45;

	QLabel *lb_ol_c_52;
	QLineEdit *le_ol_c_52;
	QLabel *lb_ol_phi_52;
	QLineEdit *le_ol_phi_52;
	QLabel *lb_ol_c_54;
	QLineEdit *le_ol_c_54;
	QLabel *lb_ol_phi_54;
	QLineEdit *le_ol_phi_54;
	QLabel *lb_ol_c_56;
	QLineEdit *le_ol_c_56;
	QLabel *lb_ol_phi_56;
	QLineEdit *le_ol_phi_56;

	MGroupBox_H *gbh_ol_ar;
	QLabel *lb_ol_ar_min;
	QLineEdit *le_ol_ar_min;
	QLabel *lb_ol_ar_max;
	QLineEdit *le_ol_ar_max;

	QWidget *wg_cond_lens;
	QWidget *wg_obj_lens;
	QTabWidget *tb_lens;

	QDockWidget *dw_microscope;

	/*******************************************************/
	mt::System_Configuration system_conf;

	mt::Atom_Data<double> atoms;

	mt::Vector<complex<double>, mt::e_host> iw_psi;
	mt::Grid_2d<double> grid_2d_iw_psi;

	mt::Detector<double, mt::e_host> stem_detector;

	mt::Atomic_Data atomic_data;

	void read_system_configuration();

	mt::Output_Multislice<float> output_multislice_float;
	mt::Output_Multislice<double> output_multislice_double;

	template<class T>
	void read_input_multislice(mt::Input_Multislice<T> &input_multislice);

	/*******************************************************/
	void set_enable_running_simulation(bool enable);

	void set_objective_lens_tab_enable(boolean enable);

	/*******************************************************/
	void default_general_dock_widget();

	void default_specimen_dock_widget();

	void default_iw_widget();
	void default_stem_widget();
	void default_pcs_widget();
	void default_eels_widget();
	void default_pptf_widget();
	void default_experiment_dock_widget();

	void default_condenser_lens_tab();
	void default_objective_lens_tab();
	void default_microscope_dock_widget();

	/*******************************************************/
	void create_iw_widget();
	void create_stem_widget();
	void create_pcs_widget();
	void create_eels_widget();
	void create_pptf_widget();
	void create_condenser_lens_tab(QWidget *&lens);
	void create_objective_lens_tab(QWidget *&lens);

	/*******************************************************/
	void create_general_dock_widget();
	void create_specimen_dock_widget();
	void create_experiment_dock_widget();
	void create_microscope_dock_widget();

	/*******************************************************/
	void set_specimen_enable(bool enable);

	/*******************************************************/
	void disable_item_cb_iw_type(mt::eIncident_Wave_Type iwt, bool disabled=true);
	void disable_item_cb_illu_model(mt::eIllumination_Model ill_model, bool disabled=true);

	/*******************************************************/
	void set_cb_iw_type_using_eTEM_Sim_Type(mt::eTEM_Sim_Type sim_type);
	void set_cb_illu_model_using_eTEM_Sim_Type(mt::eTEM_Sim_Type sim_type);

};

/************************************************/
class RS_Thread : public QThread
{
	Q_OBJECT

public:
	RS_Thread(MainWindow *main_i): main(main_i){}

signals:
	void resultReady();

protected:
	void run()
	{
		main->run_sim();
		emit resultReady();
	}
private:
	MainWindow *main;
};

class PB_Timer : public QObject
{
 Q_OBJECT

public:
	PB_Timer(MainWindow *main_i): main(main_i)
	{
		connect(&timer, SIGNAL (timeout()), this, SLOT (doWork()));
		timer.setInterval(500);
	}

	void start()
	{
		timer.start();
	}

	void stop()
	{
		timer.stop();
	}
private slots:
	void doWork()
	{
		main->progress_run_update();
	}

private:
	MainWindow *main;
	QTimer timer;
};
#endif // MAINWINDOW_H
