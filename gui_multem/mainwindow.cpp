#include "types.cuh"
#include "memory_info.cuh"
#include "lin_alg_def.cuh"
#include "atomic_data.hpp"
#include "atomic_data_mt.hpp"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "tem_simulation.cuh"
#include "multem.cu"

#include <QtWidgets>
#include <QThread>

#include "q_types.h"
#include "q_widget.h"
#include "q_load_specimen.h"
#include "q_data_viewer.h"
#include "mainwindow.h"

MainWindow::MainWindow(): central_widget(new QLabel)
{
	rs_thread = new RS_Thread(this);
	connect(rs_thread, SIGNAL(resultReady()), this, SLOT(show_data()));

	pb_timer = new PB_Timer(this);

	// set seed to generate cgpu_rand numbers
	qsrand(QTime::currentTime().msec());

	setWindowTitle(tr("MULTEM"));
	setDockOptions(QMainWindow::AnimatedDocks);

	create_menubar();

	central_widget->setFrameStyle(QFrame::Box | QFrame::Sunken);
	central_widget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::MinimumExpanding);

	setCentralWidget(central_widget);

	create_general_dock_widget();
	create_specimen_dock_widget();
	create_experiment_dock_widget();
	create_microscope_dock_widget();

	statusBar()->showMessage("MULTEM");

	QList<QComboBox*> allComboBoxes = findChildren<QComboBox *>();
	for (auto icb = 0; icb < allComboBoxes.size(); icb++)
	{
		connect(allComboBoxes[icb], static_cast<void(QComboBox::*)(int)>(&QComboBox::highlighted),
				[=](int index){statusBar()->showMessage(allComboBoxes[icb]->itemData(index, Qt::ToolTipRole).toString());});
	}

	QList<QWidget*> allQWidgets= findChildren<QWidget *>();
	for (auto iw = 0; iw < allQWidgets.size(); iw++)
	{
		allQWidgets[iw]->setToolTip(allQWidgets[iw]->statusTip());
	}

	connect(cb_iw_type, SIGNAL(currentIndexChanged(int)), SLOT(cb_iw_type_currentIndexChanged(int)));

	connect(cb_sim_type, SIGNAL(currentIndexChanged(int)), SLOT(cb_sim_type_currentIndexChanged(int)));
	auto sim_type = cb_sim_type->findData(mt::eTEMST_EWRS);
	cb_sim_type->setCurrentIndex(sim_type);
	cb_sim_type_currentIndexChanged(sim_type);

	connect(cb_elec_spec_int_model, SIGNAL(currentIndexChanged(int)), SLOT(cb_elec_spec_int_model_currentIndexChanged(int)));
	cb_elec_spec_int_model->setCurrentIndex(0);
	cb_elec_spec_int_model_currentIndexChanged(0);

	set_specimen_enable(false);
	default_sim();

	//  connect(cb_potential_type, static_cast<void(QComboBox::*)(int)>(&QComboBox::highlighted),
	//	[&](int index){statusBar()->showMessage(cb_potential_type->itemData(index, Qt::ToolTipRole).toString());});

	//  QList<QComboBox*> allPButtons = findChildren<QComboBox *>();
	//  qDebug() << "Amount of children found :" << allPButtons.count();

	//	[=](int index){statusBar()->showMessage(itemData(index, Qt::ToolTipRole).toString());};
	//  auto cc = QColor("#F1D8D8").toRgb();
	//  setStyleSheet("QFrame { background: "+ cc.name() +" }");
	//  setStyleSheet("QPushButton, QLineEdit, QComboBox { color: red }");
	//  setStyleSheet("QComboBox::drop-down { image: url(dropdown.png) }");
}

void MainWindow::open_sim()
{

}

void MainWindow::save_sim()
{

}

void MainWindow::save_as_sim()
{

}

void MainWindow::exit()
{

}

void MainWindow::about()
{

}

void MainWindow::default_sim()
{
	default_general_dock_widget();

	default_specimen_dock_widget();

	default_experiment_dock_widget();

	default_microscope_dock_widget();

	act_default_sim->setEnabled(true);
	act_start_sim->setEnabled(false);
	act_stop_sim->setEnabled(false);
}

void MainWindow::start_sim()
{
	set_enable_running_simulation(false);

	rs_thread->start();

	pb_timer->start();
}

void MainWindow::stop_sim()
{
	mt_stop_multislice(system_conf);
}

/**********************************************************/
void MainWindow::progress_run_update()
{
	auto niter = mt_niter(system_conf);
	auto iter = mt_iter(system_conf);
	pb_progress_run->setMaximum(niter);
	pb_progress_run->setValue(iter);
}

void MainWindow::run_sim()
{
	read_system_configuration();

	/**********************************************/
//	if(system_conf.is_float_host())
//	{
		mt::Input_Multislice<float> input_multislice;
		read_input_multislice(input_multislice);

		output_multislice_float.set_input_data(&input_multislice);

		mt_run_multislice<float, mt::e_host>(system_conf, input_multislice, output_multislice_float);
//	}

//	if(system_conf.is_double_host())
//	{
//		mt::Input_Multislice<double> input_multislice;
//		read_input_multislice(input_multislice);

//		mt::Output_Multislice<double> output_multislice;
//		output_multislice.set_input_data(&input_multislice);

//		run_multislice<double, mt::e_host>(system_conf, input_multislice, output_multislice);

//		if(success_multislice<double, mt::e_host>())
//		{
//			auto image = new Data_Viewer(this);
//			image->set_input_data(output_multislice);
//			image->show();
//		}
//	}

//	if(system_conf.is_float_device())
//	{
//		mt::Input_Multislice<float> input_multislice;
//		read_input_multislice(input_multislice);

//		mt::Output_Multislice<float> output_multislice;
//		output_multislice.set_input_data(&input_multislice);

//		run_multislice<float, mt::e_device>(system_conf, input_multislice, output_multislice);

//		if(success_multislice<float, mt::e_device>())
//		{
//			auto image = new Data_Viewer(this);
//			image->set_input_data(output_multislice);
//			image->show();
//		}
//	}

//	if(system_conf.is_double_device())
//	{
//		mt::Input_Multislice<double> input_multislice;
//		read_input_multislice(input_multislice);

//		mt::Output_Multislice<double> output_multislice;
//		output_multislice.set_input_data(&input_multislice);

//		run_multislice<double, mt::e_device>(system_conf, input_multislice, output_multislice);

//		if(success_multislice<double, mt::e_device>())
//		{
//			auto image = new Data_Viewer(this);
//			image->set_input_data(output_multislice);
//			image->show();
//		}
//	}
}

void MainWindow::show_data()
{
	pb_timer->stop();

	pb_progress_run->setValue(pb_progress_run->maximum());

	set_enable_running_simulation(true);

	if(mt_success_multislice(system_conf))
	{
		statusBar()->showMessage("Simulation ready");

		auto image = new Data_Viewer(output_multislice_float);
		image->show();
		image->activateWindow();
	}
	else
	{
		output_multislice_float.clear();
		pb_progress_run->setMinimum(0);
		pb_progress_run->setMaximum(100);
		pb_progress_run->setValue(0);
	}
}
/**********************************************************/
void MainWindow::cb_device_currentIndexChanged(int index)
{
	auto device = cb_device->itemData(index);
	auto bb_cpu = (device == mt::e_host);

	lb_nthreads->setVisible(bb_cpu);
	sb_nthreads->setVisible(bb_cpu);
	lb_gpu_card->setVisible(!bb_cpu);
	cb_gpu_card->setVisible(!bb_cpu);
}

void MainWindow::cb_elec_spec_int_model_currentIndexChanged(int index)
{
	auto elec_spec_int_model = cb_elec_spec_int_model->itemData(index);
	auto bb = (elec_spec_int_model == mt::eESIM_Multislice);
	gbg_potential_slic->setEnabled(bb);
}

void MainWindow::pb_spec_load_released()
{
	QString filename = QFileDialog::getOpenFileName(this,
		tr("Open specimen txt file"), "", tr("Specimen (*.txt)"));

	QFileInfo info(filename);
	default_specimen_dock_widget();

	if(!filename.isEmpty() && info.isFile())
	{
		Load_Specimen<double> load_specimen;
		bool success = false;

		try
		{
			success = load_specimen(filename, atoms);
			statusBar()->showMessage("Specimen loaded", 1000);
		}
		catch(...)
		{
			atoms.clear();
			success = false;
			le_spec_file->clear();
			lb_pptf_n_slices->setText(tr(n_slices_qba(0)));
			cb_pptf_slice->clear();
			wg_pptf->setEnabled(false);

			statusBar()->showMessage("An error occurred while trying to load the specimen file", 1000);
		}

		set_specimen_enable(success);
		if(success)
		{
			le_spec_file->setText(filename);

			auto n_atoms = atoms.size();
			auto n_types = atoms.Z_unique.size();

			lb_spec_n_atom_types->setText(tr(n_atoms_types_qba(n_atoms, n_types)));
			le_spec_lx->setText(QString::number(atoms.l_x, 'f', 4));
			le_spec_ly->setText(QString::number(atoms.l_y, 'f', 4));
			pb_spec_recenter->setEnabled(false);

			le_thk_0->setText(QString::number(atoms.z_min, 'f', 4));
			le_thk_d->setText(QString::number(atoms.dz, 'f', 4));
			le_thk_e->setText(QString::number(atoms.z_max, 'f', 4));

			le_pot_slic_thick->setText(QString::number(atoms.dz, 'f', 4));

			/****************************************************/
			le_iw_x->setText(QString::number(le_spec_lx->text().toDouble()/2, 'f', 3));
			le_iw_y->setText(QString::number(le_spec_ly->text().toDouble()/2, 'f', 3));

			/****************************************************/
			le_stem_sc_px_e->setText(QString::number(le_spec_lx->text().toDouble(), 'f', 3));
			le_stem_sc_py_e->setText(QString::number(le_spec_ly->text().toDouble(), 'f', 3));

			disconnect(cb_eels_element, SIGNAL(currentIndexChanged(int)), this, SLOT(cb_eels_element_currentIndexChanged(int)));

			cb_eels_element->clear();
			for(auto itype=0; itype<n_types; itype++)
			{
				auto Z = atoms.Z_unique[itype];
				cb_eels_element->addItem(QString::fromStdString(atomic_data.Z_name(Z)), QVariant(Z));
			}
			cb_eels_element_currentIndexChanged(0);

			connect(cb_eels_element, SIGNAL(currentIndexChanged(int)), SLOT(cb_eels_element_currentIndexChanged(int)));

		}
		act_start_sim->setEnabled(success);
	}
	else
	{
		atoms.clear();
		set_specimen_enable(false);
		le_spec_file->clear();
		le_stem_sc_px_e->setText("0.000");
		le_stem_sc_py_e->setText("0.000");
		cb_eels_element->clear();

		lb_pptf_n_slices->setText(tr(n_slices_qba(0)));
		cb_pptf_slice->clear();
		wg_pptf->setEnabled(false);

		act_start_sim->setEnabled(false);
	}
}

void MainWindow::le_spec_lx_ly_textChanged(const QString &)
{
	pb_spec_recenter->setEnabled(atoms.size()>0);
	le_iw_x->setText(QString::number(le_spec_lx->text().toDouble()/2, 'f', 3));
	le_iw_y->setText(QString::number(le_spec_ly->text().toDouble()/2, 'f', 3));

}

void MainWindow::pb_spec_recenter_released()
{
	if(atoms.size()>0)
	{
		auto lx = le_spec_lx->text().toDouble();
		auto ly = le_spec_ly->text().toDouble();

		atoms.xy_recenter(lx, ly);

		statusBar()->showMessage("The specimen was recenter along x-y directions");
	}
}

void MainWindow::cb_spec_rot_center_type_currentIndexChanged(int index)
{
	auto spec_rot_p0_type = cb_spec_rot_center_type->itemData(index);
	auto bb = !(spec_rot_p0_type == mt::eRPT_geometric_center);
	lb_spec_rot_center_p->setEnabled(bb);
	le_spec_rot_center_px->setEnabled(bb);
	le_spec_rot_center_py->setEnabled(bb);
	le_spec_rot_center_pz->setEnabled(bb);
}

void MainWindow::cb_pot_slic_type_currentIndexChanged(int index)
{
	auto potential_slic = cb_pot_slic_type->itemData(index);
	auto bb = (potential_slic != mt::ePS_Planes)&&(potential_slic != mt::ePS_Auto);
	lb_pot_slic_thick->setEnabled(bb);
	le_pot_slic_thick->setEnabled(bb);
}

void MainWindow::cb_thk_type_currentIndexChanged(int index)
{
	auto thk_type = cb_thk_type->itemData(index);
	auto bb_ws = (thk_type==mt::eTT_Whole_Spec);
	lb_thk_0->setVisible(!bb_ws);
	le_thk_0->setVisible(!bb_ws);
	lb_thk_d->setVisible(!bb_ws);
	le_thk_d->setVisible(!bb_ws);
	lb_thk_e->setVisible(!bb_ws);
	le_thk_e->setVisible(!bb_ws);
}

void MainWindow::le_pot_slic_thick_editingFinished()
{
	atoms.dz = le_pot_slic_thick->text().toDouble();
}

void MainWindow::cb_sim_type_currentIndexChanged(int index)
{
	auto sim_type = QVar_to_enum<mt::eTEM_Sim_Type>(cb_sim_type->itemData(index));

	wg_electron_phonon->setEnabled(!mt::is_IWFS_IWRS(sim_type));

	/***********************************************************/
	set_cb_iw_type_using_eTEM_Sim_Type(sim_type);
	auto iw_type = QVar_to_enum<mt::eIncident_Wave_Type>(cb_iw_type->currentData());
	iw_type = mt::validate_incident_wave_type(sim_type, iw_type);

	/***********************************************************/
	auto bb_ss = (!mt::is_specimen_required(sim_type))?true:false;
	if(bb_ss)
	{
		le_spec_file->clear();
		default_specimen_dock_widget();
		set_specimen_enable(false);

		le_spec_lx->setText(QString::number(40, 'f', 4));
		le_spec_ly->setText(QString::number(40, 'f', 4));

		le_iw_x->setText(QString::number(le_spec_lx->text().toDouble()/2, 'f', 3));
		le_iw_y->setText(QString::number(le_spec_ly->text().toDouble()/2, 'f', 3));

		pb_spec_load->setEnabled(false);
	}
	else
	{
		bb_ss = !le_spec_file->text().isEmpty();
		pb_spec_load->setEnabled(true);
	}
	act_start_sim->setEnabled(bb_ss);

	lb_spec_lx->setEnabled(bb_ss);
	le_spec_lx->setEnabled(bb_ss);
	lb_spec_ly->setEnabled(bb_ss);
	le_spec_ly->setEnabled(bb_ss);

	/***********************************************************/
	auto bb_cl = mt::is_scanning(sim_type)||mt::is_CBED_CBEI(sim_type)||mt::is_IWFS_IWRS(sim_type);
	bb_cl = bb_cl||mt::is_EWFS_EWRS(sim_type)||mt::is_EFTEM(sim_type);
	bb_cl = bb_cl && mt::is_convergent_wave(iw_type);

	wg_cond_lens->setEnabled(bb_cl);

	/***********************************************************/
	auto bb_ob = mt::is_ISTEM(sim_type)||mt::is_CBEI(sim_type);
	bb_ob = bb_ob||mt::is_HRTEM(sim_type)||mt::is_HCTEM(sim_type)||mt::is_EFTEM(sim_type);

	set_objective_lens_tab_enable(bb_ob);

	auto bb_ob_defocus  =bb_ob||mt::is_EWFS_EWRS(sim_type);

	lb_ol_sa_c_10_zero->setEnabled(bb_ob_defocus);
	cb_ol_sa_c_10_zero->setEnabled(bb_ob_defocus);

	bb_ob_defocus = bb_ob_defocus && (cb_ol_sa_c_10_zero->currentData()==mt::eZDT_User_Define);

	lb_ol_sa_c_10_z->setEnabled(bb_ob_defocus);
	le_ol_sa_c_10_z->setEnabled(bb_ob_defocus);

	/***********************************************************/
	wg_stem->setVisible(false);
	wg_pcs->setVisible(false);
	wg_eels->setVisible(false);
	wg_pptf->setVisible(false);

	if(mt::is_STEM_ISTEM(sim_type))
	{
		wg_stem->setVisible(true);
		gbg_stem_detector->setVisible(mt::is_STEM(sim_type));
	}
	else if(mt::is_PED_HCTEM(sim_type))
	{
		auto str = (mt::is_PED(sim_type))?tr("Precession"):tr("Hollow cone");
		gbg_pcs->setTitle(str);
		wg_pcs->setVisible(true);
	}
	else if(mt::is_EELS_EFTEM(sim_type))
	{
		wg_eels->setVisible(true);
		auto str = (mt::is_EELS(sim_type))?tr("Collection angle"):tr("Objective aperture");
		le_eels_coll_angle->setStatusTip(str);
	}
	else if(mt::is_PPFS_PPRS(sim_type) || mt::is_TFFS_TFRS(sim_type))
	{
		auto str = (mt::is_PPFS_PPRS(sim_type))?tr("Projected potential"):tr("Transmission function");
		gbh_pptf->setTitle(str);
		wg_pptf->setEnabled(!le_spec_file->text().isEmpty());
		wg_pptf->setVisible(true);
	}

	// select obj lens
	auto idx_lens = (mt::is_plane_wave(iw_type))?1:0;
	tb_lens->setCurrentIndex(idx_lens);

	// hide incident wave
	auto bb_iw = !(mt::is_PPFS_PPRS(sim_type)||mt::is_TFFS_TFRS(sim_type));
	wg_iw->setVisible(bb_iw);

	set_cb_illu_model_using_eTEM_Sim_Type(sim_type);
}

void MainWindow::cb_iw_type_currentIndexChanged(int index)
{
	auto sim_type = QVar_to_enum<mt::eTEM_Sim_Type>(cb_sim_type->currentData());
	auto iw_type = QVar_to_enum<mt::eIncident_Wave_Type>(cb_iw_type->itemData(index));
	iw_type = mt::validate_incident_wave_type(sim_type, iw_type);

	/***********************************************************/
	wg_cond_lens->setEnabled(mt::is_convergent_wave(iw_type));

	/***********************************************************/
	auto bb_ud = mt::is_user_define_wave(iw_type);
	lb_iw_file->setVisible(bb_ud);
	le_iw_file->setVisible(bb_ud);
	pb_iw_load->setVisible(bb_ud);

	auto bb_pw = (!mt::is_scanning(sim_type) && mt::is_convergent_wave(iw_type)) || bb_ud;

	lb_iw_x->setVisible(bb_pw);
	le_iw_x->setVisible(bb_pw);
	lb_iw_y->setVisible(bb_pw);
	le_iw_y->setVisible(bb_pw);
	pb_iw_p->setVisible(bb_pw);

	if(!bb_ud)
	{
		iw_psi.clear();
		grid_2d_iw_psi.set_input_data(0, 0, 0, 0);
		le_iw_file->clear();
	}

	// select obj lens
	auto idx_lens = (mt::is_convergent_wave(iw_type))?0:1;
	tb_lens->setCurrentIndex(idx_lens);
}

void MainWindow::pb_iw_load_released()
{
	QString filename = QFileDialog::getOpenFileName(this,
		tr("Open user incident wave binary file"), "", tr("Incident wave (*.bin)"));

	QFileInfo info(filename);

	if(!filename.isEmpty() && info.isFile())
	{
		try
		{
			mt::read_mat_binary_matrix(filename.toUtf8(), grid_2d_iw_psi, iw_psi);

			QString str = "User define incident wave loaded (nx, ny) = (" + QString::number(grid_2d_iw_psi.nx);
			str += ", " + QString::number(grid_2d_iw_psi.ny) + ")";
			statusBar()->showMessage(str);

			le_iw_file->setText(filename);
		}
		catch(...)
		{
			le_iw_file->clear();
			statusBar()->showMessage("An error occurred while trying to load the user define incident wave");
		}
	}
	else
	{
		le_iw_file->clear();
	}
}

void MainWindow::cb_stem_det_type_currentIndexChanged(int index)
{
	stem_detector.type = QVar_to_enum<mt::eDetector_Type>(cb_stem_det_type->itemData(index));
	auto bb_det_type = stem_detector.is_detector_circular();

	fr_stem_detector_ideal->setVisible(bb_det_type);
	fr_stem_detector_exp->setVisible(!bb_det_type);

	stem_detector.clear();
	stem_detector.resize(le_stem_det_n_det->text().toInt());

	cb_stem_det_k_det->setCurrentIndex(0);

	le_stem_det_ideal_inner_angle->setText("40.000");
	le_stem_det_ideal_outer_angle->setText("100.000");

	le_stem_det_exp_file->clear();

	if(bb_det_type)
	{
		stem_detector.g_inner[0] = le_stem_det_ideal_inner_angle->text().toDouble();
		stem_detector.g_outer[0] = le_stem_det_ideal_outer_angle->text().toDouble();
	}
}

void MainWindow::le_stem_det_n_det_editingFinished()
{
	auto n_det = le_stem_det_n_det->text().toInt();

	if(n_det==stem_detector.size())
	{
		return;
	}

	stem_detector.type = QVar_to_enum<mt::eDetector_Type>(cb_stem_det_type->currentData());
	stem_detector.clear();
	stem_detector.resize(n_det);

	le_stem_det_ideal_inner_angle->setText("40.000");
	le_stem_det_ideal_outer_angle->setText("100.000");

	auto bb_det_type = stem_detector.is_detector_circular();
	if(bb_det_type)
	{
		stem_detector.g_inner[0] = le_stem_det_ideal_inner_angle->text().toDouble();
		stem_detector.g_outer[0] = le_stem_det_ideal_outer_angle->text().toDouble();
	}

	disconnect(cb_stem_det_k_det, SIGNAL(currentIndexChanged(int)), this, SLOT(cb_stem_det_k_det_currentIndexChanged(int)));

	cb_stem_det_k_det->clear();
	for(auto i_det = 0; i_det<n_det; i_det++)
	{
		cb_stem_det_k_det->addItem(QString::number(i_det+1), QVariant(i_det));
	}
	cb_stem_det_k_det->setCurrentIndex(0);

	connect(cb_stem_det_k_det, SIGNAL(currentIndexChanged(int)), SLOT(cb_stem_det_k_det_currentIndexChanged(int)));

}

void MainWindow::cb_stem_det_k_det_currentIndexChanged(int index)
{
	switch (stem_detector.type)
	{
		case mt::eDT_Circular:
		{
			auto inner_ang = stem_detector.g_inner[index];
			le_stem_det_ideal_inner_angle->setText(QString::number(inner_ang, 'f', 3));
			auto outer_ang = stem_detector.g_outer[index];
			le_stem_det_ideal_outer_angle->setText(QString::number(outer_ang, 'f', 3));
		}
		break;
		case mt::eDT_Radial:
		{
//			auto fx = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "fx");
//			mt::assign(fx, input_multislice.detector.fx[i]);
		}
		break;
		case mt::eDT_Matrix:
		{
			le_stem_det_exp_file->setText(QString::fromStdString(stem_detector.fn[index]));
		}
		break;
	}
}

void MainWindow::le_stem_det_ideal_editingFinished()
{
	auto idx = cb_stem_det_k_det->currentData().toInt();
	stem_detector.g_inner[idx] = le_stem_det_ideal_inner_angle->text().toDouble();
	stem_detector.g_outer[idx] = le_stem_det_ideal_outer_angle->text().toDouble();
}

void MainWindow::pb_stem_det_exp_load_released()
{
	QString filename = QFileDialog::getOpenFileName(this,
		tr("Open detector binary file"), "", tr("Detector file (*.bin)"));

	QFileInfo info(filename);

	if(!filename.isEmpty() && info.isFile())
	{
		try
		{
			auto idx = cb_stem_det_k_det->currentData().toInt();
			stem_detector.fn[idx] = filename.toStdString();
			mt::read_mat_binary_matrix(filename.toUtf8(), stem_detector.grid_2d[idx], stem_detector.fR[idx]);

			QString str = "Detector sensitivity loaded (nx, ny) = (" + QString::number(stem_detector.grid_2d[idx].nx);
			str += ", " + QString::number(stem_detector.grid_2d[idx].ny) + ")";
			statusBar()->showMessage(str);

			le_stem_det_exp_file->setText(filename);
		}
		catch(...)
		{
			le_stem_det_exp_file->clear();
			statusBar()->showMessage("An error occurred while trying to load the user define incident wave");
		}
	}
	else
	{
		le_stem_det_exp_file->clear();
	}
}

void MainWindow::cb_eels_element_currentIndexChanged(int index)
{
	int Z = cb_eels_element->itemData(index).toInt();

	cb_eels_energy->clear();
	if(Z>0)
	{
		auto energy = atomic_data.extract_major_edges(Z);
		for(auto ik=0; ik<energy.size(); ik++)
		{
			cb_eels_energy->addItem(QString::number(energy[ik]));
		}
	}
}

void MainWindow::cb_illu_model_currentIndexChanged(int index)
{
	auto im_model = cb_illu_model->itemData(index);
	auto bb_incoh = (im_model != mt::eIM_Coherent);
	auto bb_fint = (im_model == mt::eIM_Full_Integration);

	lb_illu_beta->setVisible(bb_incoh);
	le_illu_beta->setVisible(bb_incoh);
	le_illu_beta->setMaximumWidth((bb_fint)?QWIDGETSIZE_MAX:80);

	lb_illu_nbeta->setVisible(bb_fint);
	le_illu_nbeta->setVisible(bb_fint);

	/***************************************************/
	lb_cl_sa_delta->setVisible(bb_incoh);
	le_cl_sa_delta->setVisible(bb_incoh);

	lb_cl_sa_ndelta->setVisible(bb_fint);
	cb_cl_sa_ndelta->setVisible(bb_fint);

	/***************************************************/
	lb_ol_sa_delta->setVisible(bb_incoh);
	le_ol_sa_delta->setVisible(bb_incoh);

	lb_ol_sa_ndelta->setVisible(bb_fint);
	cb_ol_sa_ndelta->setVisible(bb_fint);
}

void MainWindow::cb_cl_sa_c_10_zero_currentIndexChanged(int index)
{
	auto c_10_z = cb_cl_sa_c_10_zero->itemData(index);
	auto bb_z = (c_10_z != mt::eZDT_First);
	lb_cl_sa_c_10_z->setEnabled(bb_z);
	le_cl_sa_c_10_z->setEnabled(bb_z);
}

void MainWindow::cb_ol_sa_c_10_zero_currentIndexChanged(int index)
{
	auto c_10_z = cb_ol_sa_c_10_zero->itemData(index);
	auto bb_z = (c_10_z == mt::eZDT_User_Define);
	lb_ol_sa_c_10_z->setEnabled(bb_z);
	le_ol_sa_c_10_z->setEnabled(bb_z);
}

void MainWindow::pb_cl_sa_c_10_opt_released()
{
	auto E0 = le_mic_acc_vol->text().toDouble();
	auto c_30 = le_cl_sa_c_30->text().toDouble()*mt::c_mm_2_Angs;
	auto f = mt::get_Scherzer_defocus(E0, c_30);
	le_cl_sa_c_10->setText(QString::number(f, 'f', 3));
}

void MainWindow::pb_ol_sa_c_10_opt_released()
{
	auto E0 = le_mic_acc_vol->text().toDouble();
	auto c_30 = le_ol_sa_c_30->text().toDouble()*mt::c_mm_2_Angs;
	auto f = mt::get_Scherzer_defocus(E0, c_30);
	le_ol_sa_c_10->setText(QString::number(f, 'f', 3));
}

void MainWindow::disable_item_cb_iw_type(mt::eIncident_Wave_Type iwt, bool disabled)
{
	auto index_dis = cb_iw_type->findData(QVariant(iwt));
	if(disabled && (index_dis==cb_iw_type->currentIndex()))
	{
		auto index_sel = cb_iw_type->findData(QVariant(mt::eIWT_Auto));
		cb_iw_type->setCurrentIndex(index_sel);
	}
	disable_item_QComboBox(cb_iw_type, index_dis, disabled);
}

void MainWindow::disable_item_cb_illu_model(mt::eIllumination_Model ill_model, bool disabled)
{
	auto index_dis = cb_illu_model->findData(QVariant(ill_model));
	if(disabled && (index_dis==cb_illu_model->currentIndex()))
	{
		auto index_sel = cb_illu_model->findData(QVariant(mt::eIM_Coherent));
		cb_illu_model->setCurrentIndex(index_sel);
	}
	disable_item_QComboBox(cb_illu_model, index_dis, disabled);
}

void MainWindow::create_menubar()
{
	file_menu = menuBar()->addMenu(tr("&File"));
	QToolBar *file_tool_bar = addToolBar(tr("File"));

	const QIcon icon_open_sim = QIcon::fromTheme("document-open", QIcon(":/images/open.png"));
	act_open_sim = new QAction(icon_open_sim, tr("&Open simulation"), this);
	act_open_sim->setShortcuts(QKeySequence::Open);
	act_open_sim->setStatusTip(tr("Open a TEM simulation"));
	connect(act_open_sim, &QAction::triggered, this, &MainWindow::open_sim);
	file_menu->addAction(act_open_sim);
	file_tool_bar->addAction(act_open_sim);
	act_open_sim->setEnabled(false);

	file_menu->addSeparator();

	const QIcon icon_save_sim = QIcon::fromTheme("document-open", QIcon(":/images/save.png"));
	act_save_sim = new QAction(icon_save_sim, tr("&Save simulation"), this);
	act_save_sim->setShortcuts(QKeySequence::Save);
	act_save_sim->setStatusTip(tr("Save the current TEM simulation"));
	connect(act_save_sim, &QAction::triggered, this, &MainWindow::save_sim);
	file_menu->addAction(act_save_sim);
	file_tool_bar->addAction(act_save_sim);
	act_save_sim->setEnabled(false);

	file_tool_bar->addSeparator();

	/******************************************************************/
	act_save_as_sim = new QAction(tr("&Save simulation as"), this);
	act_save_as_sim->setShortcuts(QKeySequence::SaveAs);
	act_save_as_sim->setStatusTip(tr("Save the current TEM simulation"));
	connect(act_save_as_sim, &QAction::triggered, this, &MainWindow::save_as_sim);
	file_menu->addAction(act_save_as_sim);

	file_menu->addSeparator();

	act_exit = new QAction(tr("&Exit"), this);
	act_exit->setShortcuts(QKeySequence::Quit);
	act_exit->setStatusTip(tr("Exit the application"));
	connect(act_exit, &QAction::triggered, this, &MainWindow::exit);
	file_menu->addAction(act_exit);

	/******************************************************************/
	const QIcon icon_default_sim = QIcon::fromTheme("Set default values", QIcon(":/images/default.png"));
	act_default_sim = new QAction(icon_default_sim, tr("&Default values"), this);
	act_default_sim->setShortcuts(QKeySequence::Save);
	act_default_sim->setStatusTip(tr("Set default values"));
	connect(act_default_sim, &QAction::triggered, this, &MainWindow::default_sim);
	file_tool_bar->addAction(act_default_sim);

	const QIcon icon_start_sim = QIcon::fromTheme("Start simulation", QIcon(":/images/start.png"));
	act_start_sim = new QAction(icon_start_sim, tr("&Start Simulation"), this);
	act_start_sim->setStatusTip(tr("Start simulation"));
	act_start_sim->setEnabled(false);
	connect(act_start_sim, &QAction::triggered, this, &MainWindow::start_sim);
	file_tool_bar->addAction(act_start_sim);

	const QIcon icon_stop_sim = QIcon::fromTheme("Stop simulation", QIcon(":/images/stop.png"));
	act_stop_sim = new QAction(icon_stop_sim, tr("&Stop Simulation"), this);
	act_stop_sim->setStatusTip(tr("Stop simulation"));
	act_stop_sim->setEnabled(false);
	connect(act_stop_sim, &QAction::triggered, this, &MainWindow::stop_sim);
	file_tool_bar->addAction(act_stop_sim);

	file_tool_bar->addSeparator();

	pb_progress_run = new QProgressBar;
	pb_progress_run->setMinimum(0);
	pb_progress_run->setMaximum(100);
	pb_progress_run->setValue(0);

	file_tool_bar->addWidget(pb_progress_run);
}

/*******************************************************/
void MainWindow::read_system_configuration()
{
	system_conf.device = QVar_to_enum<mt::eDevice>(cb_device->currentData());
	system_conf.precision = QVar_to_enum<mt::ePrecision>(cb_precision->currentData());
	system_conf.cpu_ncores = 1;
	system_conf.cpu_nthread = sb_nthreads->value();
	system_conf.gpu_device = cb_gpu_card->currentData().toInt();
	system_conf.gpu_nstream = 1;
	system_conf.active = true;

	system_conf.validate_parameters();
}

template<class T>
void MainWindow::read_input_multislice(mt::Input_Multislice<T> &input_multislice)
{
	/************************ simulation type **************************/
	input_multislice.simulation_type = QVar_to_enum<mt::eTEM_Sim_Type>(cb_sim_type->currentData());

	/*******************************************************************/
	input_multislice.interaction_model = QVar_to_enum<mt::eElec_Spec_Int_Model>(cb_elec_spec_int_model->currentData());
	input_multislice.potential_type = QVar_to_enum<mt::ePotential_Type>(cb_potential_type->currentData());

	/*******************************************************************/
	input_multislice.operation_mode = mt::eOM_Normal;
	input_multislice.reverse_multislice = false;

	/************** Electron-Phonon interaction model ******************/
	input_multislice.pn_model = QVar_to_enum<mt::ePhonon_Model>(cb_pn_model->currentData());
	input_multislice.pn_coh_contrib = ckb_pn_coh_contrib->isChecked();
	input_multislice.pn_single_conf = ckb_pn_single_conf->isChecked();
	input_multislice.pn_nconf = le_pn_nconf->text().toInt();
	input_multislice.pn_dim.set(cb_pn_dim->currentData().toInt());
	input_multislice.pn_seed = le_pn_seed->text().toInt();

	if(input_multislice.is_specimen_required())
	{
		/**************************** Specimen *****************************/
		auto ct_na = 1;
		auto ct_nb = 1;
		auto ct_nc = 1;
		auto ct_a = 0;
		auto ct_b = 0;
		auto ct_c = 0;
		auto ct_x0 = 0;
		auto ct_y0 = 0;

		mt::Vector<mt::Amorp_Lay_Info<T>, mt::e_host> amorp_lay_info(0);
		for(auto i = 0; i<amorp_lay_info.size(); i++)
		{
			amorp_lay_info[i].z_0 = 0;
			amorp_lay_info[i].z_e = 0;
			amorp_lay_info[i].dz = 2.0;
		}

		input_multislice.atoms.set_crystal_parameters(ct_na, ct_nb, ct_nc, ct_a, ct_b, ct_c, ct_x0, ct_y0);
		input_multislice.atoms.set_amorphous_parameters(amorp_lay_info);
		input_multislice.atoms.set_atoms(atoms);

		/************************ Specimen rotation *************************/
		input_multislice.spec_rot_theta = le_spec_rot_theta->text().toDouble()*mt::c_deg_2_rad;
		auto u_x = le_spec_rot_ux->text().toDouble();
		auto u_y = le_spec_rot_uy->text().toDouble();
		auto u_z = le_spec_rot_uz->text().toDouble();
		input_multislice.spec_rot_u0 = mt::r3d<T>(u_x, u_y, u_z);
		input_multislice.spec_rot_center_type = QVar_to_enum<mt::eRot_Point_Type>(cb_spec_rot_center_type->currentData());
		auto p_x = le_spec_rot_center_px->text().toDouble();
		auto p_y = le_spec_rot_center_py->text().toDouble();
		auto p_z = le_spec_rot_center_pz->text().toDouble();
		input_multislice.spec_rot_center_p = mt::r3d<T>(p_x, p_y, p_z);

		/************************ Specimen thickness ***********************/
		input_multislice.thick_type = QVar_to_enum<mt::eThick_Type>(cb_thk_type->currentData());
		if(!input_multislice.is_whole_spec())
		{
			auto thk_0 = le_thk_0->text().toDouble();
			auto thk_d = le_thk_d->text().toDouble();
			auto thk_e = le_thk_e->text().toDouble() + mt::epsilon_rel<float>();

			auto n_thk = int(floor((thk_e-thk_0)/thk_d+0.5))+2;
			input_multislice.thick.clear();
			input_multislice.thick.reserve(n_thk);

			auto thk = thk_0;
			while(thk<=thk_e)
			{
				input_multislice.thick.push_back(thk);
				thk += thk_d;
			}
		}

		/************************ Potential slicing ************************/
		input_multislice.potential_slicing = QVar_to_enum<mt::ePotential_Slicing>(cb_pot_slic_type->currentData());
//		atoms.dz = le_pot_slic_thick->text().toDouble();
//		input_multislice.atoms.dz = atoms.dz;
	}

	auto lx = le_spec_lx->text().toDouble();
	auto ly = le_spec_ly->text().toDouble();
	auto dz = le_pot_slic_thick->text().toDouble();

	/************************** xy sampling ****************************/
	auto nx = le_samp_nx->text().toDouble();
	auto ny = le_samp_nx->text().toDouble();
	bool bwl = ckb_samp_bandwidth->isChecked();
	bool pbc_xy = true;

	input_multislice.grid_2d.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	/************************ Incident wave ****************************/
	auto iw_type = QVar_to_enum<mt::eIncident_Wave_Type>(cb_iw_type->currentData());
	input_multislice.set_incident_wave_type(iw_type);

	if(input_multislice.is_user_define_wave())
	{
		mt::assign(iw_psi, input_multislice.iw_psi);
	}

	// read iw_x and iw_y
	auto iw_x = le_iw_x->text().toDouble();
	auto iw_y = le_iw_y->text().toDouble();

	input_multislice.iw_x.clear();
	input_multislice.iw_x.push_back(iw_x);

	input_multislice.iw_y.clear();
	input_multislice.iw_y.push_back(iw_y);

	/********************* Microscope parameter ***********************/
	input_multislice.E_0 = le_mic_acc_vol->text().toDouble();
	input_multislice.theta = le_tilt_theta->text().toDouble();
	input_multislice.phi = le_tilt_phi->text().toDouble();

	/********************* Illumination model *************************/
	input_multislice.illumination_model = QVar_to_enum<mt::eIllumination_Model>(cb_illu_model->currentData());
	input_multislice.temporal_spatial_incoh = mt::eTSI_Temporal_Spatial;

	input_multislice.cond_lens.beta = le_illu_beta->text().toDouble()*mt::c_mrad_2_rad; 			// divergence semi-angle (mrad-->rad)
	input_multislice.cond_lens.nbeta = le_illu_nbeta->text().toInt();								// Number of integration steps for the divergence semi-angle

	/*********************** Objective lens **************************/
	/********************* Symmetric aberrations *********************/
	input_multislice.cond_lens.c_30 = le_cl_sa_c_30->text().toDouble()*mt::c_mm_2_Angs; 			// 3rd order spherical aberration (mm-->Angstrom)
	input_multislice.cond_lens.c_50 = le_cl_sa_c_50->text().toDouble()*mt::c_mm_2_Angs;			 // 5th order spherical aberration (mm-->Angstrom)
	input_multislice.cond_lens.m = le_cl_sa_m->text().toInt();									  // momentum of the vortex

	input_multislice.cond_lens.c_10 = le_cl_sa_c_10->text().toDouble();	 							// defocus (Angstrom)
	input_multislice.cond_lens.zero_defocus_type = QVar_to_enum<mt::eZero_Defocus_Type>(cb_cl_sa_c_10_zero->currentData());	 // Zero defocus type
	input_multislice.cond_lens.zero_defocus_plane = le_cl_sa_c_10_z->text().toDouble();										 // Zero defocus position
	input_multislice.cond_lens.sf = le_cl_sa_delta->text().toDouble(); 								// Defocus spread (Angstrom)
	input_multislice.cond_lens.nsf = cb_cl_sa_ndelta->currentData().toInt(); 								// Number of integration steps for the defocus Spread

	/******************** Non-symmetric aberrations ******************/
	input_multislice.cond_lens.c_12 = le_cl_c_12->text().toDouble();	 							// 2-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_12 = le_cl_phi_12->text().toDouble()*mt::c_deg_2_rad;	 		// Azimuthal angle of 2-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.c_21 = le_cl_c_21->text().toDouble(); 								// Axial coma (Angstrom)
	input_multislice.cond_lens.phi_21 = le_cl_phi_21->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of axial coma (degrees-->rad)
	input_multislice.cond_lens.c_23 = le_cl_c_23->text().toDouble(); 								// 3-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_23 = le_cl_phi_23->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 3-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.c_32 = le_cl_c_32->text().toDouble(); 								// Axial star aberration (Angstrom)
	input_multislice.cond_lens.phi_32 = le_cl_phi_32->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of axial star aberration (degrees-->rad)
	input_multislice.cond_lens.c_34 = le_cl_c_34->text().toDouble(); 								// 4-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_34 = le_cl_phi_34->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 4-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.c_41 = le_cl_c_41->text().toDouble();  								// 4th order axial coma (Angstrom)
	input_multislice.cond_lens.phi_41 = le_cl_phi_41->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 4th order axial coma (degrees-->rad)
	input_multislice.cond_lens.c_43 = le_cl_c_43->text().toDouble(); 								// 3-lobe aberration (Angstrom)
	input_multislice.cond_lens.phi_43 = le_cl_phi_43->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 3-lobe aberration (degrees-->rad)
	input_multislice.cond_lens.c_45 = le_cl_c_45->text().toDouble(); 								// 5-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_45 = le_cl_phi_45->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 5-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.c_52 = le_cl_c_52->text().toDouble();								// 5th order axial star aberration (Angstrom)
	input_multislice.cond_lens.phi_52 = le_cl_phi_52->text().toDouble()*mt::c_deg_2_rad;			// Azimuthal angle of 5th order axial star aberration (degrees-->rad)
	input_multislice.cond_lens.c_54 = le_cl_c_54->text().toDouble();								// 5th order rosette aberration (Angstrom)
	input_multislice.cond_lens.phi_54 = le_cl_phi_54->text().toDouble()*mt::c_deg_2_rad;			// Azimuthal angle of 5th order rosette aberration (degrees-->rad)
	input_multislice.cond_lens.c_56 = le_cl_c_56->text().toDouble();								// 6-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_56 = le_cl_phi_56->text().toDouble()*mt::c_deg_2_rad;			// Azimuthal angle of 6-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.inner_aper_ang = le_cl_ar_min->text().toDouble()*mt::c_mrad_2_rad;	// inner aperture (mrad-->rad)
	input_multislice.cond_lens.outer_aper_ang = le_cl_ar_max->text().toDouble()*mt::c_mrad_2_rad;	// outer aperture (mrad-->rad)

	input_multislice.cond_lens.set_input_data(input_multislice.E_0, input_multislice.grid_2d);

	/*********************** Objective lens **************************/
	/********************* Symmetric aberrations *********************/
	input_multislice.obj_lens.c_30 = le_ol_sa_c_30->text().toDouble()*mt::c_mm_2_Angs;			  // 3rd order spherical aberration (mm-->Angstrom)
	input_multislice.obj_lens.c_50 = le_ol_sa_c_50->text().toDouble()*mt::c_mm_2_Angs;			 // 5th order spherical aberration (mm-->Angstrom)
	input_multislice.obj_lens.m = le_ol_sa_m->text().toInt();									  // momentum of the vortex

	input_multislice.obj_lens.c_10 = le_ol_sa_c_10->text().toDouble();	 							// defocus (Angstrom)
	input_multislice.obj_lens.zero_defocus_type = QVar_to_enum<mt::eZero_Defocus_Type>(cb_ol_sa_c_10_zero->currentData());	 // Zero defocus type
	input_multislice.obj_lens.zero_defocus_plane = le_ol_sa_c_10_z->text().toDouble();										 // Zero defocus position
	input_multislice.obj_lens.sf = le_ol_sa_delta->text().toDouble(); 								// Defocus spread (Angstrom)
	input_multislice.obj_lens.nsf = cb_ol_sa_ndelta->currentData().toInt(); 								// Number of integration steps for the defocus Spread

	/****************** Non-symmetric aberrations *******************/
	input_multislice.obj_lens.c_12 = le_ol_c_12->text().toDouble();								 // 2-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_12 = le_ol_phi_12->text().toDouble()*mt::c_deg_2_rad;			 // Azimuthal angle of 2-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.c_21 = le_ol_c_21->text().toDouble(); 								// Axial coma (Angstrom)
	input_multislice.obj_lens.phi_21 = le_ol_phi_21->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of axial coma (degrees-->rad)
	input_multislice.obj_lens.c_23 = le_ol_c_23->text().toDouble(); 								// 3-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_23 = le_ol_phi_23->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 3-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.c_32 = le_ol_c_32->text().toDouble(); 								// Axial star aberration (Angstrom)
	input_multislice.obj_lens.phi_32 = le_ol_phi_32->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of axial star aberration (degrees-->rad)
	input_multislice.obj_lens.c_34 = le_ol_c_34->text().toDouble(); 								// 4-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_34 = le_ol_phi_34->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 4-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.c_41 = le_ol_c_41->text().toDouble();  								// 4th order axial coma (Angstrom)
	input_multislice.obj_lens.phi_41 = le_ol_phi_41->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 4th order axial coma (degrees-->rad)
	input_multislice.obj_lens.c_43 = le_ol_c_43->text().toDouble(); 								// 3-lobe aberration (Angstrom)
	input_multislice.obj_lens.phi_43 = le_ol_phi_43->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 3-lobe aberration (degrees-->rad)
	input_multislice.obj_lens.c_45 = le_ol_c_45->text().toDouble(); 								// 5-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_45 = le_ol_phi_45->text().toDouble()*mt::c_deg_2_rad; 			// Azimuthal angle of 5-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.c_52 = le_ol_c_52->text().toDouble();								 // 5th order axial star aberration (Angstrom)
	input_multislice.obj_lens.phi_52 = le_ol_phi_52->text().toDouble()*mt::c_deg_2_rad;			 // Azimuthal angle of 5th order axial star aberration (degrees-->rad)
	input_multislice.obj_lens.c_54 = le_ol_c_54->text().toDouble();								 // 5th order rosette aberration (Angstrom)
	input_multislice.obj_lens.phi_54 = le_ol_phi_54->text().toDouble()*mt::c_deg_2_rad;			 // Azimuthal angle of 5th order rosette aberration (degrees-->rad)
	input_multislice.obj_lens.c_56 = le_ol_c_56->text().toDouble();								 // 6-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_56 = le_ol_phi_56->text().toDouble()*mt::c_deg_2_rad;			 // Azimuthal angle of 6-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.inner_aper_ang = le_ol_ar_min->text().toDouble()*mt::c_mrad_2_rad;	// inner aperture (mrad-->rad)
	input_multislice.obj_lens.outer_aper_ang = le_ol_ar_max->text().toDouble()*mt::c_mrad_2_rad;	// outer aperture (mrad-->rad)

	input_multislice.obj_lens.set_input_data(input_multislice.E_0, input_multislice.grid_2d);

	/************************* ISTEM/STEM **************************/
	if(input_multislice.is_scanning())
	{
		input_multislice.scanning.type = QVar_to_enum<mt::eScanning_Type>(cb_stem_sc_type->currentData());
		input_multislice.scanning.pbc = ckb_stem_sc_incl_last_point->isChecked();
		input_multislice.scanning.ns = le_stem_sc_n_points->text().toInt();
		input_multislice.scanning.x0 = le_stem_sc_px_0->text().toDouble();
		input_multislice.scanning.y0 = le_stem_sc_py_0->text().toDouble();
		input_multislice.scanning.xe = le_stem_sc_px_e->text().toDouble();
		input_multislice.scanning.ye = le_stem_sc_py_e->text().toDouble();
		input_multislice.scanning.set_grid();
	}

	if(input_multislice.is_STEM())
	{
		input_multislice.detector.type = QVar_to_enum<mt::eDetector_Type>(cb_stem_det_type->currentData());
		input_multislice.detector.assign(stem_detector);

		if(input_multislice.detector.is_detector_circular())
		{
			auto lambda = mt::get_lambda(input_multislice.E_0);
			for(auto idet=0; idet<input_multislice.detector.size(); idet++)
			{
				auto g_inner = input_multislice.detector.g_inner[idet];
				input_multislice.detector.g_inner[idet] = sin(g_inner*mt::c_mrad_2_rad)/lambda;

				auto g_outer = input_multislice.detector.g_outer[idet];
				input_multislice.detector.g_outer[idet] = sin(g_outer*mt::c_mrad_2_rad)/lambda;
			}
		}
	}
	else if (input_multislice.is_PED_HCTEM())
	{
		input_multislice.theta = le_pcs_tilt_angle->text().toDouble()*mt::c_deg_2_rad;
		input_multislice.nrot = le_pcs_n_rotation->text().toInt();
	}
	else if (input_multislice.is_EELS_EFTEM())
	{
		mt::eSpace space = (input_multislice.is_EELS())?mt::eS_Reciprocal:mt::eS_Real;

		auto Z = cb_eels_element->currentData().toInt();
		auto E_loss = cb_eels_energy->currentText().toDouble()*mt::c_eV_2_keV;
		auto collection_angle = le_eels_coll_angle->text().toDouble()*mt::c_mrad_2_rad;
		auto m_selection = cb_eels_m_selection->currentData().toInt();
		auto channelling_type = QVar_to_enum<mt::eChannelling_Type>(cb_eels_m_selection->currentData());

		input_multislice.eels_fr.set_input_data(space, input_multislice.E_0, E_loss, m_selection, collection_angle, channelling_type, Z);
	}

	input_multislice.validate_parameters();
}

void MainWindow::set_specimen_enable(bool enable)
{
	lb_spec_n_atom_types->setEnabled(enable);
	lb_spec_lx->setEnabled(enable);
	le_spec_lx->setEnabled(enable);
	lb_spec_ly->setEnabled(enable);
	le_spec_ly->setEnabled(enable);
	pb_spec_recenter->setEnabled(enable);

	gbg_spec_rot->setEnabled(enable);
	gbg_thk->setEnabled(enable);
	gbg_potential_slic->setEnabled(enable);
}

void MainWindow::set_cb_iw_type_using_eTEM_Sim_Type(mt::eTEM_Sim_Type sim_type)
{
	auto iw_type = (mt::is_IWFS_IWRS(sim_type))?mt::eIWT_Convergent_Wave:mt::eIWT_Auto;
	auto idx_iw_type = cb_iw_type->findData(QVariant(iw_type));
	cb_iw_type->setCurrentIndex(idx_iw_type);
	cb_iw_type_currentIndexChanged(idx_iw_type);

	// disable plane wave
	auto bb_plane = mt::is_scanning(sim_type)||mt::is_CBED_CBEI(sim_type);
	disable_item_cb_iw_type(mt::eIWT_Plane_Wave, bb_plane);

	// disable convergent beam
	auto bb_cb = mt::is_ED_HRTEM(sim_type)||mt::is_PED_HCTEM(sim_type);
	disable_item_cb_iw_type(mt::eIWT_Convergent_Wave, bb_cb);

	// disable user define wave
	auto bb_user_def = mt::is_PED_HCTEM(sim_type);
	disable_item_cb_iw_type(mt::eIWT_User_Define_Wave, bb_user_def);

	// disable auto wave
	auto bb_auto = mt::is_IWFS_IWRS(sim_type);
	disable_item_cb_iw_type(mt::eIWT_Auto, bb_auto);
}

void MainWindow::set_cb_illu_model_using_eTEM_Sim_Type(mt::eTEM_Sim_Type sim_type)
{
	auto bb_full_int = mt::is_EELS_EFTEM(sim_type)||mt::is_EWFS_EWRS(sim_type);
	bb_full_int = bb_full_int||mt::is_IWFS_IWRS(sim_type)||mt::is_PPFS_PPRS(sim_type);
	bb_full_int = bb_full_int||mt::is_TFFS_TFRS(sim_type);

	auto bb_part_coh = bb_full_int||mt::is_STEM(sim_type)||mt::is_CBED(sim_type);

	disable_item_cb_illu_model(mt::eIM_Coherent, false);
	disable_item_cb_illu_model(mt::eIM_Partial_Coherent, bb_part_coh);
	disable_item_cb_illu_model(mt::eIM_Trans_Cross_Coef, true);
	disable_item_cb_illu_model(mt::eIM_Full_Integration, bb_full_int);

	cb_illu_model->setCurrentIndex(cb_illu_model->findData(mt::eIM_Coherent));
}

/*******************************************************/
void MainWindow::default_general_dock_widget()
{
	if(cb_gpu_card->count()>0)
	{
		cb_device->setCurrentIndex(cb_device->findData(mt::e_device));
		cb_gpu_card->setCurrentIndex(0);
	}
	else
	{
		cb_device->setCurrentIndex(cb_device->findData(mt::e_host));
	}

	cb_precision->setCurrentIndex(cb_precision->findData(mt::eP_float));

	sb_nthreads->setValue(sb_nthreads->maximum());

	cb_elec_spec_int_model->setCurrentIndex(cb_elec_spec_int_model->findData(mt::eESIM_Multislice));
	cb_potential_type->setCurrentIndex(cb_potential_type->findData(mt::ePT_Lobato_0_12));

	cb_pn_model->setCurrentIndex(cb_pn_model->findData(mt::ePM_Frozen_Phonon));
	ckb_pn_coh_contrib->setChecked(false);
	ckb_pn_single_conf->setChecked(false);
	le_pn_nconf->setText("5");
	cb_pn_dim->setCurrentIndex(cb_pn_dim->findData(110));
	le_pn_seed->setText("300183");
	ckb_pn_auto_seed->setChecked(false);
}

void MainWindow::default_specimen_dock_widget()
{
	/*******************************************************/
	/**********************Specimen info********************/
	/*******************************************************/
	le_spec_file->clear();

	lb_spec_n_atom_types->setText(tr(n_atoms_types_qba(0, 0)));
	le_spec_lx->setText("0.000");
	le_spec_ly->setText("0.000");
	pb_spec_recenter->setEnabled(false);

	/*******************************************************/
	/*********************** Rotation **********************/
	/*******************************************************/
	le_spec_rot_theta->setText("0.000");
	le_spec_rot_ux->setText("1.0");
	le_spec_rot_uy->setText("0.0");
	le_spec_rot_uz->setText("0.0");

	/*******************************************************/
	cb_spec_rot_center_type->setCurrentIndex(0);
	cb_spec_rot_center_type_currentIndexChanged(0);

	le_spec_rot_center_px->setText("0.0");
	le_spec_rot_center_py->setText("0.0");
	le_spec_rot_center_pz->setText("0.0");

	/*******************************************************/
	/*********************** Thickness *********************/
	/*******************************************************/
	cb_thk_type->setCurrentIndex(0);
	cb_thk_type_currentIndexChanged(0);

	/***************************************************/
	le_thk_0->setText("0.0000");
	le_thk_d->setText("0.0000");
	le_thk_e->setText("0.0000");

	/*******************************************************/
	/***************** Potential slicing *******************/
	/*******************************************************/
	cb_pot_slic_type->setCurrentIndex(0);
	cb_pot_slic_type_currentIndexChanged(0);

	le_pot_slic_thick->setText("2.0000");

	/*******************************************************/
	/********************** xy sampling ********************/
	/*******************************************************/
	le_samp_nx->setText("1024");
	le_samp_ny->setText("1024");
	ckb_samp_bandwidth->setChecked(false);

	/*******************************************************/
	/******************** Clear atoms **********************/
	/*******************************************************/
	atoms.clear();

	le_iw_x->setText("0.000");
	le_iw_y->setText("0.000");
	le_stem_sc_px_e->setText("0.000");
	le_stem_sc_py_e->setText("0.000");

	cb_eels_energy->clear();
	cb_eels_element->clear();
}

void MainWindow::default_iw_widget()
{
	cb_iw_type->setCurrentIndex(cb_iw_type->findData(mt::eIWT_Auto));
	le_iw_x->setText("0.000");
	le_iw_y->setText("0.000");
	le_iw_file->clear();
}

void MainWindow::default_stem_widget()
{
	cb_stem_sc_type->setCurrentIndex(1);
	ckb_stem_sc_incl_last_point->setChecked(true);
	le_stem_sc_n_points->setText("10");
	cb_stem_sc_np_assign->setCurrentIndex(0);
	le_stem_sc_px_0->setText("0.000");
	le_stem_sc_py_0->setText("0.000");
	le_stem_sc_px_e->setText("0.000");
	le_stem_sc_py_e->setText("0.000");

	cb_stem_det_type->setCurrentIndex(cb_stem_det_type->findData(mt::eDT_Circular));
	le_stem_det_n_det->setText("1");
	cb_stem_det_type->setCurrentIndex(0);
	cb_stem_det_type_currentIndexChanged(0);

	stem_detector.resize(cb_stem_det_k_det->count());
	le_stem_det_ideal_inner_angle->setText("40.000");
	le_stem_det_ideal_outer_angle->setText("100.000");

	le_stem_det_exp_file->clear();
}

void MainWindow::default_pcs_widget()
{
	le_pcs_tilt_angle->setText("1.5");
	le_pcs_n_rotation->setText("4");

}

void MainWindow::default_eels_widget()
{
	cb_eels_energy->clear();
	cb_eels_element->clear();
	le_eels_coll_angle->setText("25.000");
	cb_eels_m_selection->setCurrentIndex(0);
	cb_eels_channelling_type->setCurrentIndex(0);
}

void MainWindow::default_pptf_widget()
{
	cb_pptf_slice->clear();
}

void MainWindow::default_experiment_dock_widget()
{
	default_iw_widget();
	default_stem_widget();
	default_pcs_widget();
	default_eels_widget();
	default_pptf_widget();

	cb_sim_type->setCurrentIndex(cb_sim_type->findData(mt::eTEMST_EWRS));
}

void MainWindow::default_condenser_lens_tab()
{
	le_cl_sa_c_30->setText("0.001");
	le_cl_sa_c_50->setText("0.000");
	le_cl_sa_m->setText("0");
	le_cl_sa_c_10->setText("0.000");
	cb_cl_sa_c_10_zero->setCurrentIndex(cb_cl_sa_c_10_zero->findData(mt::eZDT_First));
	le_cl_sa_c_10_z->setText("0.000");
	le_cl_sa_delta->setText("20.0");
	cb_cl_sa_ndelta->setCurrentIndex(6);
	pb_cl_sa_c_10_opt_released();

	le_cl_c_12->setText("0.000");
	le_cl_phi_12->setText("0.000");

	le_cl_c_21->setText("0.000");
	le_cl_phi_21->setText("0.000");

	le_cl_c_23->setText("0.000");
	le_cl_phi_23->setText("0.000");

	le_cl_c_32->setText("0.000");
	le_cl_phi_32->setText("0.000");

	le_cl_c_34->setText("0.000");
	le_cl_phi_34->setText("0.000");

	le_cl_c_41->setText("0.000");
	le_cl_phi_41->setText("0.000");

	le_cl_c_43->setText("0.000");
	le_cl_phi_43->setText("0.000");

	le_cl_c_45->setText("0.000");
	le_cl_phi_45->setText("0.000");

	le_cl_c_52->setText("0.000");
	le_cl_phi_52->setText("0.000");

	le_cl_c_54->setText("0.000");
	le_cl_phi_54->setText("0.000");

	le_cl_c_56->setText("0.000");
	le_cl_phi_56->setText("0.000");

	le_cl_ar_min->setText("0.000");
	le_cl_ar_max->setText("21.000");
}

void MainWindow::default_objective_lens_tab()
{
	le_ol_sa_c_30->setText("0.001");
	le_ol_sa_c_50->setText("0.000");
	le_ol_sa_m->setText("0");
	le_ol_sa_c_10->setText("0.000");
	cb_ol_sa_c_10_zero->setCurrentIndex(cb_ol_sa_c_10_zero->findData(mt::eZDT_Last));
	le_ol_sa_c_10_z->setText("0.000");
	le_ol_sa_delta->setText("20.0");
	cb_ol_sa_ndelta->setCurrentIndex(6);
	pb_ol_sa_c_10_opt_released();

	le_ol_c_12->setText("0.000");
	le_ol_phi_12->setText("0.000");

	le_ol_c_21->setText("0.000");
	le_ol_phi_21->setText("0.000");

	le_ol_c_23->setText("0.000");
	le_ol_phi_23->setText("0.000");

	le_ol_c_32->setText("0.000");
	le_ol_phi_32->setText("0.000");

	le_ol_c_34->setText("0.000");
	le_ol_phi_34->setText("0.000");

	le_ol_c_41->setText("0.000");
	le_ol_phi_41->setText("0.000");

	le_ol_c_43->setText("0.000");
	le_ol_phi_43->setText("0.000");

	le_ol_c_45->setText("0.000");
	le_ol_phi_45->setText("0.000");

	le_ol_c_52->setText("0.000");
	le_ol_phi_52->setText("0.000");

	le_ol_c_54->setText("0.000");
	le_ol_phi_54->setText("0.000");

	le_ol_c_56->setText("0.000");
	le_ol_phi_56->setText("0.000");

	le_ol_ar_min->setText("0.000");
	le_ol_ar_max->setText("0.000");
}

void MainWindow::default_microscope_dock_widget()
{
	cb_mic_name->setCurrentIndex(0);
	le_mic_acc_vol->setText("300.0");

	le_tilt_theta->setText("0.000");
	le_tilt_phi->setText("0.000");

	cb_illu_model->setCurrentIndex(cb_illu_model->findData(mt::eIM_Coherent));
	le_illu_beta->setText("0.100");
	le_illu_nbeta->setText("8");

	default_condenser_lens_tab();
	default_objective_lens_tab();
}

/*******************************************************/
void MainWindow::set_enable_running_simulation(bool enable)
{
	act_default_sim->setEnabled(enable);
	act_start_sim->setEnabled(enable);
	act_stop_sim->setEnabled(!enable);

	dw_general->setEnabled(enable);
	dw_specimen->setEnabled(enable);
	dw_simulation->setEnabled(enable);
	dw_microscope->setEnabled(enable);
}

void MainWindow::set_objective_lens_tab_enable(boolean enable)
{
	lb_ol_sa_c_30->setEnabled(enable);
	le_ol_sa_c_30->setEnabled(enable);
	lb_ol_sa_c_50->setEnabled(enable);
	le_ol_sa_c_50->setEnabled(enable);
	lb_ol_sa_m->setEnabled(enable);
	le_ol_sa_m->setEnabled(enable);

	lb_ol_sa_c_10->setEnabled(enable);
	le_ol_sa_c_10->setEnabled(enable);
	pb_ol_sa_c_10_opt->setEnabled(enable);
	lb_ol_sa_c_10_zero->setEnabled(enable);
	cb_ol_sa_c_10_zero->setEnabled(enable);
	lb_ol_sa_c_10_z->setEnabled(enable);
	le_ol_sa_c_10_z->setEnabled(enable);

	lb_ol_sa_delta->setEnabled(enable);
	le_ol_sa_delta->setEnabled(enable);
	lb_ol_sa_ndelta->setEnabled(enable);
	cb_ol_sa_ndelta->setEnabled(enable);

	gbge_ol_nsa->setEnabled(enable);

	gbh_ol_ar->setEnabled(enable);
}

/*******************************************************/
void MainWindow::create_iw_widget()
{
	lb_iw_type = new QLabel(tr("Type"));

	cb_iw_type = new QComboBox;
	cb_iw_type->setStatusTip(tr("Incident wave type"));
	cb_iw_type->addItem(tr("Plane wave"), QVariant(mt::eIWT_Plane_Wave));
	cb_iw_type->setItemData(0, tr("Plane wave illumination"), Qt::ToolTipRole);
	cb_iw_type->addItem(tr("Convergent wave"), QVariant(mt::eIWT_Convergent_Wave));
	cb_iw_type->setItemData(1, tr("Convergent wave illumination"), Qt::ToolTipRole);
	cb_iw_type->addItem(tr("User define"), QVariant(mt::eIWT_User_Define_Wave));
	cb_iw_type->setItemData(2, tr("User define (matrix)"), Qt::ToolTipRole);
	cb_iw_type->addItem(tr("Auto"), QVariant(mt::eIWT_Auto));
	cb_iw_type->setItemData(3, tr("It is defined by the experiment"), Qt::ToolTipRole);

	/***************************************************/
	lb_iw_x = new QLabel(tr("X<sub>0</sub> [Å]"));

	le_iw_x = new QLineEdit;
	le_iw_x->setStatusTip(tr("Initial x scan position"));
	le_iw_x->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_iw_x->setText("0.000");

	lb_iw_y = new QLabel(tr("Y<sub>0</sub> [Å]"));

	le_iw_y = new QLineEdit;
	le_iw_y->setStatusTip(tr("Initial y scan position"));
	le_iw_y->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_iw_y->setText("0.000");

	pb_iw_p = new QPushButton(QIcon(":/images/cross_pointer.png"), "");
	pb_iw_p->setStatusTip(tr("Select position of the incidente wave"));
	pb_iw_p->setMaximumSize(22, 22);

	/***************************************************/
	lb_iw_file = new QLabel(tr("Filename"));

	le_iw_file = new QLineEdit;
	le_iw_file->setStatusTip(tr("User define incident wave filename"));
	le_iw_file->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);

	set_non_editable_QLineEdit(le_iw_file);

	pb_iw_load = new QPushButton(QIcon(":/images/open.png"), "");
	pb_iw_load->setStatusTip(tr("Load user define incident wave"));
	pb_iw_load->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	/***************************************************/
	gbg_iw = new MGroupBox_G(tr("Incident wave"));
	gbg_iw->addWidget(lb_iw_type, 0, 0);
	gbg_iw->addWidget(cb_iw_type, 0, 1, 1, 4);

	gbg_iw->addWidget(lb_iw_x, 1, 0);
	gbg_iw->addWidget(le_iw_x, 1, 1);
	gbg_iw->addWidget(lb_iw_y, 1, 2);
	gbg_iw->addWidget(le_iw_y, 1, 3);
	gbg_iw->addWidget(pb_iw_p, 1, 4);

	gbg_iw->addWidget(lb_iw_file, 2, 0);
	gbg_iw->addWidget(le_iw_file, 2, 1, 1, 3);
	gbg_iw->addWidget(pb_iw_load, 2, 4);

	connect(pb_iw_load, SIGNAL(released()), SLOT(pb_iw_load_released()));

	auto layout = new QVBoxLayout;
	layout->setSpacing(0);
	layout->setContentsMargins(0, 0, 0, 0);
	layout->setSizeConstraint(QLayout::SetMinimumSize);
	layout->addWidget(gbg_iw);

	wg_iw = new QWidget;
	wg_iw->setLayout(layout);
}

void MainWindow::create_stem_widget()
{
	/*******************************************************/
	/**********************scannning************************/
	/*******************************************************/

	lb_stem_sc_type = new QLabel(tr("Type"));

	cb_stem_sc_type = new QComboBox;
	cb_stem_sc_type->setStatusTip(tr("Scan type"));
	cb_stem_sc_type->addItem(tr("Line"), QVariant(1));
	cb_stem_sc_type->setItemData(0, tr("Scan along a line"), Qt::ToolTipRole);
	cb_stem_sc_type->addItem(tr("Area"), QVariant(2));
	cb_stem_sc_type->setItemData(1, tr("Scan an area"), Qt::ToolTipRole);
	cb_stem_sc_type->setCurrentIndex(1);

	lb_stem_sc_n_points = new QLabel(tr(x_sub_to_qba("N", "scans", "x-large")));
	lb_stem_sc_n_points->setAlignment(Qt::AlignHCenter);

	le_stem_sc_n_points = new QLineEdit;
	le_stem_sc_n_points->setMaximumWidth(45);
	le_stem_sc_n_points->setStatusTip(tr("Number of scanning points"));
	le_stem_sc_n_points->setValidator(new QIntValidator(2, int_lim));
	le_stem_sc_n_points->setText("10");

	lb_stem_sc_np_assign = new QLabel(tr("Assign"));

	cb_stem_sc_np_assign = new QComboBox;
	cb_stem_sc_np_assign->setStatusTip(tr("Assign number of scanning points to Lx or Ly"));
	cb_stem_sc_np_assign->addItem(tr("max"), QVariant(1));
	cb_stem_sc_np_assign->setItemData(0, tr("Assign number of scanning points to the maximun between Lx & Ly"), Qt::ToolTipRole);
	cb_stem_sc_np_assign->addItem(tr("Lx"), QVariant(2));
	cb_stem_sc_np_assign->setItemData(1, tr("Assign number of scanning points to Lx"), Qt::ToolTipRole);
	cb_stem_sc_np_assign->addItem(tr("Ly"), QVariant(3));
	cb_stem_sc_np_assign->setItemData(2, tr("Assign number of scanning points to Ly"), Qt::ToolTipRole);
	cb_stem_sc_np_assign->setCurrentIndex(0);

	ckb_stem_sc_incl_last_point = new QCheckBox(tr("pbc"));
	ckb_stem_sc_incl_last_point->setStatusTip(tr("Exclude last point in the scanning process"));
	ckb_stem_sc_incl_last_point->setChecked(true);

	auto lyh_scanning = new QHBoxLayout;
	lyh_scanning->setSpacing(3);
	lyh_scanning->addWidget(cb_stem_sc_type);
	lyh_scanning->addWidget(lb_stem_sc_n_points);
	lyh_scanning->addWidget(le_stem_sc_n_points);
	lyh_scanning->addWidget(cb_stem_sc_np_assign);
	lyh_scanning->addWidget(ckb_stem_sc_incl_last_point);
		/********************************************************/
	lb_stem_sc_px_0 = new QLabel(tr("X<sub>0</sub> [Å]"));

	le_stem_sc_px_0 = new QLineEdit;
	le_stem_sc_px_0->setStatusTip(tr("Initial x scan position"));
	le_stem_sc_px_0->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_stem_sc_px_0->setText("0.000");

	lb_stem_sc_py_0 = new QLabel(tr("Y<sub>0</sub> [Å]"));

	le_stem_sc_py_0 = new QLineEdit;
	le_stem_sc_py_0->setStatusTip(tr("Initial y scan position"));
	le_stem_sc_py_0->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_stem_sc_py_0->setText("0.000");

	pb_stem_sc_p_0 = new QPushButton(QIcon(":/images/cross_pointer.png"), "");
	pb_stem_sc_p_0->setStatusTip(tr("Select initial scanning point"));
	pb_stem_sc_p_0->setMaximumSize(22, 22);

	lb_stem_sc_px_e = new QLabel(tr("X<sub>e</sub> [Å]"));

	le_stem_sc_px_e = new QLineEdit;
	le_stem_sc_px_e->setStatusTip(tr("Final x scan position"));
	le_stem_sc_px_e->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_stem_sc_px_e->setText("0.000");

	lb_stem_sc_py_e = new QLabel(tr("Y<sub>e</sub> [Å]"));

	le_stem_sc_py_e = new QLineEdit;
	le_stem_sc_py_e->setStatusTip(tr("Final y scan position"));
	le_stem_sc_py_e->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_stem_sc_py_e->setText("0.000");

	pb_stem_sc_p_e = new QPushButton(QIcon(":/images/cross_pointer.png"), "");
	pb_stem_sc_p_e->setStatusTip(tr("Select end scanning point"));
	pb_stem_sc_p_e->setMaximumSize(22, 22);

	gbg_stem_scanning = new MGroupBox_G(tr("Scanning"));
	gbg_stem_scanning->layout->setVerticalSpacing(0);

	gbg_stem_scanning->layout->addLayout(lyh_scanning, 0, 0, 1, 5);

	gbg_stem_scanning->layout->setRowMinimumHeight(2, 2);

	auto hl_stem_scanning = new MLine_H;
	gbg_stem_scanning->addWidget(hl_stem_scanning, 3, 0, 1, 5);

	gbg_stem_scanning->addWidget(lb_stem_sc_px_0, 4, 0);
	gbg_stem_scanning->addWidget(le_stem_sc_px_0, 4, 1);
	gbg_stem_scanning->addWidget(lb_stem_sc_py_0, 4, 2);
	gbg_stem_scanning->addWidget(le_stem_sc_py_0, 4, 3);
	gbg_stem_scanning->addWidget(pb_stem_sc_p_0, 4, 4);

	gbg_stem_scanning->addWidget(lb_stem_sc_px_e, 5, 0);
	gbg_stem_scanning->addWidget(le_stem_sc_px_e, 5, 1);
	gbg_stem_scanning->addWidget(lb_stem_sc_py_e, 5, 2);
	gbg_stem_scanning->addWidget(le_stem_sc_py_e, 5, 3);
	gbg_stem_scanning->addWidget(pb_stem_sc_p_e, 5, 4);

	/*******************************************************/
	/***********************detectors***********************/
	/*******************************************************/
	cb_stem_det_type = new QComboBox;
	cb_stem_det_type->setStatusTip(tr("Detector types"));
	cb_stem_det_type->addItem(tr("Ideal"), mt::eDT_Circular);
	cb_stem_det_type->setItemData(0, tr("Ideal annular detector sensitivity"), Qt::ToolTipRole);
	cb_stem_det_type->addItem(tr("Radial"), mt::eDT_Radial);
	cb_stem_det_type->setItemData(1, tr("Radial detector sensitivity"), Qt::ToolTipRole);
	cb_stem_det_type->addItem(tr("Matrix"), mt::eDT_Matrix);
	cb_stem_det_type->setItemData(2, tr("Experimental image detector sensitivity"), Qt::ToolTipRole);

	disable_item_QComboBox(cb_stem_det_type, 1);

	/*******************************************************/
	lb_stem_det_n_det = new QLabel(tr(x_sub_to_qba("N", "det.", "x-large")));
	lb_stem_det_n_det->setAlignment(Qt::AlignCenter);

	le_stem_det_n_det = new QLineEdit;
	le_stem_det_n_det->setMaximumWidth(45);
	le_stem_det_n_det->setStatusTip(tr("Number of detectors"));
	le_stem_det_n_det->setValidator(new QIntValidator(1, int_lim));
	le_stem_det_n_det->setText("1");

	/*******************************************************/
	lb_stem_det_k_det = new QLabel(tr("Detector"));
	lb_stem_det_k_det->setAlignment(Qt::AlignCenter);

	cb_stem_det_k_det = new QComboBox;
	cb_stem_det_k_det->setStatusTip(tr("Select specific detector"));
	cb_stem_det_k_det->addItem(tr("1"), QVariant(0));
	cb_stem_det_k_det->setCurrentIndex(0);

	stem_detector.resize(cb_stem_det_k_det->count());
	/***************************************************/
	pb_stem_det_view = new QPushButton(QIcon(":/images/view.png"), "");
	pb_stem_det_view->setStatusTip(tr("Show Detectpr"));
	pb_stem_det_view->setMaximumSize(22, 22);
	pb_stem_det_view->setEnabled(false);

	/***************************************************/
	auto lyh_stem_det= new QHBoxLayout;
	lyh_stem_det->setContentsMargins(0, 0, 0, 0);
	lyh_stem_det->setSpacing(2);
	lyh_stem_det->addWidget(cb_stem_det_type);
	lyh_stem_det->addWidget(lb_stem_det_n_det);
	lyh_stem_det->addWidget(le_stem_det_n_det);
	lyh_stem_det->addWidget(lb_stem_det_k_det);
	lyh_stem_det->addWidget(cb_stem_det_k_det);
	lyh_stem_det->addWidget(pb_stem_det_view);

	connect(le_stem_det_n_det, SIGNAL(editingFinished()), SLOT(le_stem_det_n_det_editingFinished()));
	connect(cb_stem_det_k_det, SIGNAL(currentIndexChanged(int)), SLOT(cb_stem_det_k_det_currentIndexChanged(int)));

	/***************************************************/
	/*******************Ideal detector******************/
	/***************************************************/

	/************Idel detector - inner angle************/
	lb_stem_det_ideal_inner_angle = new QLabel(tr("Inner [mrad]"));

	le_stem_det_ideal_inner_angle = new QLineEdit;
	le_stem_det_ideal_inner_angle->setStatusTip(tr("Inner detector angle"));
	le_stem_det_ideal_inner_angle->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_stem_det_ideal_inner_angle->setText("40.000");

	/************Idel detector - outer angle************/
	lb_stem_det_ideal_outer_angle = new QLabel(tr("Outer [mrad]"));

	le_stem_det_ideal_outer_angle = new QLineEdit;
	le_stem_det_ideal_outer_angle->setStatusTip(tr("Outer detector angle"));
	le_stem_det_ideal_outer_angle->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_stem_det_ideal_outer_angle->setText("100.000");

	/***************************************************/
	auto lyh_stem_det_ideal = new QHBoxLayout;
	lyh_stem_det_ideal->setContentsMargins(2, 2, 2, 2);
	lyh_stem_det_ideal->addWidget(lb_stem_det_ideal_inner_angle);
	lyh_stem_det_ideal->addWidget(le_stem_det_ideal_inner_angle);
	lyh_stem_det_ideal->addWidget(lb_stem_det_ideal_outer_angle);
	lyh_stem_det_ideal->addWidget(le_stem_det_ideal_outer_angle);

	fr_stem_detector_ideal = new QFrame;
	fr_stem_detector_ideal->setFrameStyle(QFrame::Box | QFrame::Sunken);
	fr_stem_detector_ideal->setLineWidth(1);
	fr_stem_detector_ideal->setLayout(lyh_stem_det_ideal);

	connect(le_stem_det_ideal_inner_angle, SIGNAL(editingFinished()), SLOT(le_stem_det_ideal_editingFinished()));
	connect(le_stem_det_ideal_outer_angle, SIGNAL(editingFinished()), SLOT(le_stem_det_ideal_editingFinished()));

	/***************************************************/
	/*************Experimental detector*****************/
	/***************************************************/

	lb_stem_det_exp_file = new QLabel(tr("Filename"));

	le_stem_det_exp_file = new QLineEdit;
	le_stem_det_exp_file->setStatusTip(tr("Experimental detector sensitivity filename"));
	le_stem_det_exp_file->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);

	set_non_editable_QLineEdit(le_stem_det_exp_file);

	pb_stem_det_exp_load = new QPushButton(QIcon(":/images/open.png"), "");
	pb_stem_det_exp_load->setStatusTip(tr("Load experimental detector sensitivity"));
	pb_stem_det_exp_load->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	auto lyg_stem_detector_exp = new QGridLayout;
	lyg_stem_detector_exp->setVerticalSpacing(0);
	lyg_stem_detector_exp->setHorizontalSpacing(2);
	lyg_stem_detector_exp->setContentsMargins(2, 2, 2, 2);
	lyg_stem_detector_exp->addWidget(lb_stem_det_exp_file, 0, 0);
	lyg_stem_detector_exp->addWidget(le_stem_det_exp_file, 0, 1);
	lyg_stem_detector_exp->addWidget(pb_stem_det_exp_load, 0, 2);

	fr_stem_detector_exp = new QFrame;
	fr_stem_detector_exp->setFrameStyle(QFrame::Box | QFrame::Sunken);
	fr_stem_detector_exp->setLineWidth(1);
	fr_stem_detector_exp->setLayout(lyg_stem_detector_exp);

	connect(pb_stem_det_exp_load, SIGNAL(released()), SLOT(pb_stem_det_exp_load_released()));
	/***************************************************/
	gbg_stem_detector = new MGroupBox_G(tr("Detectors"));
	gbg_stem_detector->layout->addLayout(lyh_stem_det, 0, 0);
	gbg_stem_detector->layout->setRowMinimumHeight(1, 3);
	gbg_stem_detector->addWidget(fr_stem_detector_ideal, 2, 0);
	gbg_stem_detector->addWidget(fr_stem_detector_exp, 2, 0);

	connect(cb_stem_det_type, SIGNAL(currentIndexChanged(int)), SLOT(cb_stem_det_type_currentIndexChanged(int)));

	cb_stem_det_type->setCurrentIndex(0);
	cb_stem_det_type_currentIndexChanged(0);

	auto layout = new QVBoxLayout;
	layout->setSpacing(0);
	layout->setContentsMargins(2, 1, 1, 0);
	layout->setSizeConstraint(QLayout::SetMinimumSize);
	layout->addWidget(gbg_stem_scanning);
	layout->addWidget(gbg_stem_detector);

	wg_stem = new QWidget;
	wg_stem->setLayout(layout);
}

void MainWindow::create_pcs_widget()
{
	lb_pcs_tilt_angle = new QLabel(tr("<span style=""font-family:symbol;font-size:large"">q</span> [°]"));

	le_pcs_tilt_angle = new QLineEdit;
	le_pcs_tilt_angle->setStatusTip(tr("Precession angle"));
	le_pcs_tilt_angle->setValidator(new QDoubleValidator(-25, 25, 3));
	le_pcs_tilt_angle->setText("1.5");

	/***************************************************/
	lb_pcs_n_rotation = new QLabel(tr("# of rotations"));

	le_pcs_n_rotation = new QLineEdit;
	le_pcs_n_rotation->setStatusTip(tr("Number of rotation angles"));
	le_pcs_n_rotation->setValidator(new QIntValidator(4, int_lim));
	le_pcs_n_rotation->setText("4");

	/***************************************************/
	gbg_pcs = new MGroupBox_G(tr("Precession"));
	gbg_pcs->addWidget(lb_pcs_tilt_angle, 0, 0);
	gbg_pcs->addWidget(le_pcs_tilt_angle, 0, 1);

	gbg_pcs->addWidget(lb_pcs_n_rotation, 0, 2);
	gbg_pcs->addWidget(le_pcs_n_rotation, 0, 3);

	auto layout = new QVBoxLayout;
	layout->setSpacing(0);
	layout->setContentsMargins(2, 1, 1, 0);
	layout->setSizeConstraint(QLayout::SetMinimumSize);
	layout->addWidget(gbg_pcs);

	wg_pcs = new QWidget;
	wg_pcs->setLayout(layout);
}

void MainWindow::create_eels_widget()
{
	lb_eels_element = new QLabel(tr("Element"));
	lb_eels_element->setAlignment(Qt::AlignCenter);

	cb_eels_element = new QComboBox;
	cb_eels_element->setStatusTip(tr("Element name"));
	cb_eels_element->setMaximumWidth(55);

	/***************************************************/
	lb_eels_energy = new QLabel(tr("Energy"));

	cb_eels_energy = new QComboBox;
	cb_eels_energy->setStatusTip(tr("Energy Loss"));
	cb_eels_energy->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);

	/***************************************************/
	lb_eels_coll_angle = new QLabel(tr("<span style=""font-family:symbol;font-size:large"">q</span> [mrad]"));

	le_eels_coll_angle = new QLineEdit;
	le_eels_coll_angle->setStatusTip(tr("Collection angle"));
	le_eels_coll_angle->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_eels_coll_angle->setText("25.000");
	le_eels_coll_angle->setMaximumWidth(55);

	/***************************************************/
	lb_eels_m_selection = new QLabel(tr("m selection"));

	cb_eels_m_selection = new QComboBox;
	cb_eels_m_selection->setStatusTip(tr("MDOS selection rule"));
	cb_eels_m_selection->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
	cb_eels_m_selection->addItem(tr("Kx,Ky,Kz"), QVariant(3));
	cb_eels_m_selection->addItem(tr("Kx"), QVariant(-2));
	cb_eels_m_selection->addItem(tr("Ky"), QVariant(2));
	cb_eels_m_selection->addItem(tr("Kz"), QVariant(0));
	cb_eels_m_selection->addItem(tr("Kx-iKy"), QVariant(-1));
	cb_eels_m_selection->addItem(tr("Kx+iKy"), QVariant(1));

	cb_eels_m_selection->setCurrentIndex(0);

	/***************************************************/
	lb_eels_channelling_type = new QLabel(tr("Approx."));

	cb_eels_channelling_type = new QComboBox;
	cb_eels_channelling_type->setStatusTip(tr("Channelling approximation type"));
	cb_eels_channelling_type->addItem(tr("Single Channelling approx."), enum_to_QVar(mt::eCT_Single_Channelling));
	cb_eels_channelling_type->setItemData(0, tr("Single Channelling approximation"), Qt::ToolTipRole);
	cb_eels_channelling_type->addItem(tr("Mixed Channelling approx."), enum_to_QVar(mt::eCT_Mixed_Channelling));
	cb_eels_channelling_type->setItemData(1, tr("Mixed Channelling approximation"), Qt::ToolTipRole);
	cb_eels_channelling_type->addItem(tr("Double Channelling"), enum_to_QVar(mt::eCT_Double_Channelling));
	cb_eels_channelling_type->setItemData(2, tr("Double Channelling"), Qt::ToolTipRole);

	cb_eels_channelling_type->setCurrentIndex(0);
	/***************************************************/
	gbg_eels = new MGroupBox_G(tr("Energy Loss"));
	gbg_eels->addWidget(lb_eels_element, 0, 0);
	gbg_eels->addWidget(cb_eels_element, 0, 1);

	gbg_eels->addWidget(lb_eels_energy, 0, 2);
	gbg_eels->addWidget(cb_eels_energy, 0, 3);

	gbg_eels->addWidget(lb_eels_coll_angle, 1, 0);
	gbg_eels->addWidget(le_eels_coll_angle, 1, 1);

	gbg_eels->addWidget(lb_eels_m_selection, 1, 2);
	gbg_eels->addWidget(cb_eels_m_selection, 1, 3);

	gbg_eels->addWidget(lb_eels_channelling_type, 2, 0);
	gbg_eels->addWidget(cb_eels_channelling_type, 2, 1, 1, 3);

	auto layout = new QVBoxLayout;
	layout->setSpacing(0);
	layout->setContentsMargins(2, 1, 1, 0);
	layout->setSizeConstraint(QLayout::SetMinimumSize);
	layout->addWidget(gbg_eels);

	wg_eels = new QWidget;
	wg_eels->setLayout(layout);

	connect(cb_eels_element, SIGNAL(currentIndexChanged(int)), SLOT(cb_eels_element_currentIndexChanged(int)));

}

void MainWindow::create_pptf_widget()
{
	lb_pptf_n_slices = new QLabel(tr(n_slices_qba(0)));
	lb_pptf_n_slices->setFrameStyle(QFrame::Box | QFrame::Sunken);
	lb_pptf_n_slices->setAlignment(Qt::AlignCenter);
	lb_pptf_n_slices->setMinimumHeight(22);
	lb_pptf_n_slices->setMinimumWidth(60);

	lb_pptf_slice = new QLabel(tr("Slice"));

	cb_pptf_slice = new QComboBox;
	cb_pptf_slice->setStatusTip(tr("Slice"));
	cb_pptf_slice->setMaximumWidth(50);

	/***************************************************/
	gbh_pptf = new MGroupBox_H(tr("Projected potential"));
	gbh_pptf->addWidget(lb_pptf_n_slices);
	gbh_pptf->layout->addSpacing(4);
	gbh_pptf->addWidget(lb_pptf_slice);
	gbh_pptf->addWidget(cb_pptf_slice);
	gbh_pptf->layout->addStretch(1);

	auto layout = new QVBoxLayout;
	layout->setSpacing(0);
	layout->setContentsMargins(2, 1, 1, 0);
	layout->setSizeConstraint(QLayout::SetMinimumSize);
	layout->addWidget(gbh_pptf);

	wg_pptf = new QWidget;
	wg_pptf->setLayout(layout);
}

void MainWindow::create_condenser_lens_tab(QWidget *&lens)
{
	/*******************************************************/
	/*****************Symmetrical aberrations***************/
	/*******************************************************/

	/************3rd order spherical aberration************/
	lb_cl_sa_c_30 = new QLabel(tr(c_nm_to_qba(30, " [mm]")));
	lb_cl_sa_c_30->setAlignment(Qt::AlignCenter);

	le_cl_sa_c_30 = new QLineEdit;
	le_cl_sa_c_30->setStatusTip(tr("3rd order spherical aberration"));
	le_cl_sa_c_30->setValidator(new QDoubleValidator(-double_lim, double_lim, 5));
	le_cl_sa_c_30->setText("0.001");

	/************5th order spherical aberration************/
	lb_cl_sa_c_50 = new QLabel(tr(c_nm_to_qba(50, " [mm]")));
	lb_cl_sa_c_50->setAlignment(Qt::AlignCenter);

	le_cl_sa_c_50 = new QLineEdit;
	le_cl_sa_c_50->setStatusTip(tr("5th order spherical aberration"));
	le_cl_sa_c_50->setValidator(new QDoubleValidator(-double_lim, double_lim, 5));
	le_cl_sa_c_50->setText("0.000");

	/*******************Vortex momentum*********************/
	lb_cl_sa_m = new QLabel(tr("m"));
	lb_cl_sa_m->setAlignment(Qt::AlignCenter);

	le_cl_sa_m = new QLineEdit;
	le_cl_sa_m->setStatusTip(tr("Vortex momentum"));
	le_cl_sa_m->setValidator(new QIntValidator(-int_lim, int_lim));
	le_cl_sa_m->setText("0");

	/*************************Defocus***********************/
	lb_cl_sa_c_10 = new QLabel(tr(c_nm_to_qba(10, " [Å]")));
	lb_cl_sa_c_10->setAlignment(Qt::AlignCenter);

	le_cl_sa_c_10 = new QLineEdit;
	le_cl_sa_c_10->setStatusTip(tr("Defocus [Å]"));
	le_cl_sa_c_10->setValidator(new QDoubleValidator(-double_lim, double_lim, 3));
	le_cl_sa_c_10->setText("0.000");

	pb_cl_sa_c_10_opt = new QPushButton(QIcon(":/images/opt.png"), "");
	pb_cl_sa_c_10_opt->setStatusTip(tr("Scherzer defocus"));
	pb_cl_sa_c_10_opt->setMaximumSize(22, 22);

	/*******************************************************/
	lb_cl_sa_c_10_zero = new QLabel(tr("Focus at"));

	cb_cl_sa_c_10_zero = new QComboBox;
	cb_cl_sa_c_10_zero->setStatusTip(tr("Position to focus the incident beam"));
	cb_cl_sa_c_10_zero->addItem(tr("First atom"), QVariant(mt::eZDT_First));
	cb_cl_sa_c_10_zero->setItemData(0, tr("Focus at the first atom"), Qt::ToolTipRole);
	cb_cl_sa_c_10_zero->addItem(tr("User def."), QVariant(mt::eZDT_User_Define));
	cb_cl_sa_c_10_zero->setItemData(1, tr("Focus at user define plane"), Qt::ToolTipRole);

	/*******************************************************/
	lb_cl_sa_c_10_z = new QLabel(tr("Plane [Å]"));

	le_cl_sa_c_10_z = new QLineEdit;
	le_cl_sa_c_10_z->setStatusTip(tr("Set focus plane"));
	le_cl_sa_c_10_z->setValidator(new QDoubleValidator(-double_lim, double_lim, 3));
	le_cl_sa_c_10_z->setText("0.000");

	/*******************************************************/
	lb_cl_sa_delta = new QLabel(tr("<span style=""font-family:symbol;font-size:large"">D</span> [Å]"));
	lb_cl_sa_delta->setAlignment(Qt::AlignCenter);

	le_cl_sa_delta = new QLineEdit;
	le_cl_sa_delta->setStatusTip(tr("Defocus spread"));
	le_cl_sa_delta->setValidator(new QDoubleValidator(0, std::numeric_limits<double>::max(), 3));
	le_cl_sa_delta->setText("20.0");

	/*******************************************************/
	lb_cl_sa_ndelta = new QLabel(tr("N<sub style=""font-family:symbol;font-size:large"">D</sub>"));
	lb_cl_sa_ndelta->setAlignment(Qt::AlignCenter);

	cb_cl_sa_ndelta = new QComboBox;
	cb_cl_sa_ndelta->setStatusTip(tr("Number of integration points"));
	for(auto ik=2; ik<=128; ik++)
	{
		cb_cl_sa_ndelta->addItem(QString::number(ik), QVariant(ik));
	}
	cb_cl_sa_ndelta->setCurrentIndex(6);

	/*******************************************************/
	auto lyg_cl_sa = new QGridLayout;
	lyg_cl_sa->setHorizontalSpacing(3);
	lyg_cl_sa->setVerticalSpacing(0);
	lyg_cl_sa->setContentsMargins(2, 1, 1, 0);

	lyg_cl_sa->addWidget(lb_cl_sa_c_30, 0, 0);
	lyg_cl_sa->addWidget(le_cl_sa_c_30, 0, 1);
	lyg_cl_sa->addWidget(lb_cl_sa_c_50, 0, 2);
	lyg_cl_sa->addWidget(le_cl_sa_c_50, 0, 3);

	lyg_cl_sa->addWidget(lb_cl_sa_m, 1, 0);
	lyg_cl_sa->addWidget(le_cl_sa_m, 1, 1);
	lyg_cl_sa->addWidget(lb_cl_sa_c_10, 1, 2);

	auto lyh_c_10_opt = new QHBoxLayout;
	lyh_c_10_opt->setSpacing(0);
	lyh_c_10_opt->addWidget(le_cl_sa_c_10);
	lyh_c_10_opt->addWidget(pb_cl_sa_c_10_opt);

	lyg_cl_sa->addLayout(lyh_c_10_opt, 1, 3);

	lyg_cl_sa->addWidget(lb_cl_sa_c_10_zero, 2, 0);
	lyg_cl_sa->addWidget(cb_cl_sa_c_10_zero, 2, 1);
	lyg_cl_sa->addWidget(lb_cl_sa_c_10_z, 2, 2);
	lyg_cl_sa->addWidget(le_cl_sa_c_10_z, 2, 3);

	lyg_cl_sa->addWidget(lb_cl_sa_delta, 3, 0);
	lyg_cl_sa->addWidget(le_cl_sa_delta, 3, 1);
	lyg_cl_sa->addWidget(lb_cl_sa_ndelta, 3, 2);
	lyg_cl_sa->addWidget(cb_cl_sa_ndelta, 3, 3);

	connect(pb_cl_sa_c_10_opt, SIGNAL(released()), SLOT(pb_cl_sa_c_10_opt_released()));
	pb_cl_sa_c_10_opt->click();

	connect(cb_cl_sa_c_10_zero, SIGNAL(currentIndexChanged(int)), SLOT(cb_cl_sa_c_10_zero_currentIndexChanged(int)));

	cb_cl_sa_c_10_zero_currentIndexChanged(cb_cl_sa_c_10_zero->findData(mt::eZDT_First));

	/*******************************************************/
	/***************Nonsymmetrical aberrations**************/
	/*******************************************************/

	/*******************2-fold astigmatism******************/
	lb_cl_c_12 = new QLabel(tr(c_nm_to_qba(12, " [Å]")));

	le_cl_c_12 = new QLineEdit;
	le_cl_c_12->setStatusTip(tr("2-fold astigmatism"));
	le_cl_c_12->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_12->setText("0.000");

	lb_cl_phi_12 = new QLabel(tr(phi_nm_to_qba(12, " [°]")));

	le_cl_phi_12 = new QLineEdit;
	le_cl_phi_12->setStatusTip(tr("Azimuthal angle of 2-fold astigmatism"));
	le_cl_phi_12->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_12->setText("0.000");

	/*******************************************************/
	/***************************c2**************************/

	/***********************Axial coma**********************/
	lb_cl_c_21 = new QLabel(tr(c_nm_to_qba(21, " [Å]")));

	le_cl_c_21 = new QLineEdit;
	le_cl_c_21->setStatusTip(tr("Axial coma"));
	le_cl_c_21->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_21->setText("0.000");

	lb_cl_phi_21 = new QLabel(tr(phi_nm_to_qba(21, " [°]")));

	le_cl_phi_21 = new QLineEdit;
	le_cl_phi_21->setStatusTip(tr("Azimuthal angle of axial coma"));
	le_cl_phi_21->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_21->setText("0.000");

	/*******************3-fold astigmatism******************/
	lb_cl_c_23 = new QLabel(tr(c_nm_to_qba(23, " [Å]")));

	le_cl_c_23 = new QLineEdit;
	le_cl_c_23->setStatusTip(tr("3-fold astigmatism"));
	le_cl_c_23->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_23->setText("0.000");

	lb_cl_phi_23 = new QLabel(tr(phi_nm_to_qba(23, " [°]")));

	le_cl_phi_23 = new QLineEdit;
	le_cl_phi_23->setStatusTip(tr("Azimuthal angle of 3-fold astigmatism"));
	le_cl_phi_23->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_23->setText("0.000");

	/*******************************************************/
	/***************************c3**************************/

	/****************Axial star aberration******************/
	lb_cl_c_32 = new QLabel(tr(c_nm_to_qba(32, " [Å]")));

	le_cl_c_32 = new QLineEdit;
	le_cl_c_32->setStatusTip(tr("Axial star aberration"));
	le_cl_c_32->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_32->setText("0.000");

	lb_cl_phi_32 = new QLabel(tr(phi_nm_to_qba(32, " [°]")));

	le_cl_phi_32 = new QLineEdit;
	le_cl_phi_32->setStatusTip(tr("Azimuthal angle of axial star aberration"));
	le_cl_phi_32->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_32->setText("0.000");

	/*******************4-fold astigmatism******************/
	lb_cl_c_34 = new QLabel(tr(c_nm_to_qba(34, " [Å]")));

	le_cl_c_34 = new QLineEdit;
	le_cl_c_34->setStatusTip(tr("4-fold astigmatism"));
	le_cl_c_34->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_34->setText("0.000");

	lb_cl_phi_34 = new QLabel(tr(phi_nm_to_qba(34, " [°]")));

	le_cl_phi_34 = new QLineEdit;
	le_cl_phi_34->setStatusTip(tr("Azimuthal angle of 4-fold astigmatism"));
	le_cl_phi_34->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_34->setText("0.000");

	/*******************************************************/
	/***************************c4**************************/

	/*****************4th order axial coma******************/
	lb_cl_c_41 = new QLabel(tr(c_nm_to_qba(41, " [Å]")));

	le_cl_c_41 = new QLineEdit;
	le_cl_c_41->setStatusTip(tr("4th order axial coma"));
	le_cl_c_41->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_41->setText("0.000");

	lb_cl_phi_41 = new QLabel(tr(phi_nm_to_qba(41, " [°]")));

	le_cl_phi_41 = new QLineEdit;
	le_cl_phi_41->setStatusTip(tr("Azimuthal angle of 4th order axial coma"));
	le_cl_phi_41->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_41->setText("0.000");

	/*******************3-lobe aberration******************/
	lb_cl_c_43 = new QLabel(tr(c_nm_to_qba(43, " [Å]")));

	le_cl_c_43 = new QLineEdit;
	le_cl_c_43->setStatusTip(tr("3-lobe aberration"));
	le_cl_c_43->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_43->setText("0.000");

	lb_cl_phi_43 = new QLabel(tr(phi_nm_to_qba(43, " [°]")));

	le_cl_phi_43 = new QLineEdit;
	le_cl_phi_43->setStatusTip(tr("Azimuthal angle of 3-lobe aberration"));
	le_cl_phi_43->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_43->setText("0.000");

	/*******************5-fold astigmatism******************/
	lb_cl_c_45 = new QLabel(tr(c_nm_to_qba(45, " [Å]")));

	le_cl_c_45 = new QLineEdit;
	le_cl_c_45->setStatusTip(tr("5-fold astigmatism"));
	le_cl_c_45->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_45->setText("0.000");

	lb_cl_phi_45 = new QLabel(tr(phi_nm_to_qba(45, " [°]")));

	le_cl_phi_45 = new QLineEdit;
	le_cl_phi_45->setStatusTip(tr("Azimuthal angle of 5-fold astigmatism"));
	le_cl_phi_45->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_45->setText("0.000");

	/*******************************************************/
	/***************************c5**************************/

	/*************5th order axial star aberration***********/
	lb_cl_c_52 = new QLabel(tr(c_nm_to_qba(52, " [Å]")));

	le_cl_c_52 = new QLineEdit;
	le_cl_c_52->setStatusTip(tr("5th order axial star aberration"));
	le_cl_c_52->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_52->setText("0.000");

	lb_cl_phi_52 = new QLabel(tr(phi_nm_to_qba(52, " [°]")));

	le_cl_phi_52 = new QLineEdit;
	le_cl_phi_52->setStatusTip(tr("Azimuthal angle of 5th order axial star aberration"));
	le_cl_phi_52->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_52->setText("0.000");

	/**************5th order rosette aberration*************/
	lb_cl_c_54 = new QLabel(tr(c_nm_to_qba(54, " [Å]")));

	le_cl_c_54 = new QLineEdit;
	le_cl_c_54->setStatusTip(tr("5th order rosette aberration"));
	le_cl_c_54->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_54->setText("0.000");

	lb_cl_phi_54 = new QLabel(tr(phi_nm_to_qba(54, " [°]")));

	le_cl_phi_54 = new QLineEdit;
	le_cl_phi_54->setStatusTip(tr("Azimuthal angle of 5th order rosette aberration"));
	le_cl_phi_54->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_54->setText("0.000");

	/*******************6-fold astigmatism******************/
	lb_cl_c_56 = new QLabel(tr(c_nm_to_qba(56, " [Å]")));

	le_cl_c_56 = new QLineEdit;
	le_cl_c_56->setStatusTip(tr("6-fold astigmatism"));
	le_cl_c_56->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_c_56->setText("0.000");

	lb_cl_phi_56 = new QLabel(tr(phi_nm_to_qba(56, " [°]")));

	le_cl_phi_56 = new QLineEdit;
	le_cl_phi_56->setStatusTip(tr("Azimuthal angle of 6-fold astigmatism"));
	le_cl_phi_56->setValidator(new QDoubleValidator(-360, 360, 3));
	le_cl_phi_56->setText("0.000");

	gbge_cl_nsa = new MGroupBox_GE(tr("Nonsymmetrical aberrations"));
	gbge_cl_nsa->layout->setVerticalSpacing(0);

	gbge_cl_nsa->addWidget(lb_cl_c_12, 0, 0);
	gbge_cl_nsa->addWidget(le_cl_c_12, 0, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_12, 0, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_12, 0, 3);

	gbge_cl_nsa->addWidget(lb_cl_c_21, 1, 0);
	gbge_cl_nsa->addWidget(le_cl_c_21, 1, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_21, 1, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_21, 1, 3);
	gbge_cl_nsa->addWidget(lb_cl_c_23, 2, 0);
	gbge_cl_nsa->addWidget(le_cl_c_23, 2, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_23, 2, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_23, 2, 3);

	gbge_cl_nsa->addWidget(lb_cl_c_32, 3, 0);
	gbge_cl_nsa->addWidget(le_cl_c_32, 3, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_32, 3, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_32, 3, 3);
	gbge_cl_nsa->addWidget(lb_cl_c_34, 4, 0);
	gbge_cl_nsa->addWidget(le_cl_c_34, 4, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_34, 4, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_34, 4, 3);

	gbge_cl_nsa->addWidget(lb_cl_c_41, 5, 0);
	gbge_cl_nsa->addWidget(le_cl_c_41, 5, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_41, 5, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_41, 5, 3);
	gbge_cl_nsa->addWidget(lb_cl_c_43, 6, 0);
	gbge_cl_nsa->addWidget(le_cl_c_43, 6, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_43, 6, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_43, 6, 3);
	gbge_cl_nsa->addWidget(lb_cl_c_45, 7, 0);
	gbge_cl_nsa->addWidget(le_cl_c_45, 7, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_45, 7, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_45, 7, 3);

	gbge_cl_nsa->addWidget(lb_cl_c_52, 8, 0);
	gbge_cl_nsa->addWidget(le_cl_c_52, 8, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_52, 8, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_52, 8, 3);
	gbge_cl_nsa->addWidget(lb_cl_c_54, 9, 0);
	gbge_cl_nsa->addWidget(le_cl_c_54, 9, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_54, 9, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_54, 9, 3);
	gbge_cl_nsa->addWidget(lb_cl_c_56, 10, 0);
	gbge_cl_nsa->addWidget(le_cl_c_56, 10, 1);
	gbge_cl_nsa->addWidget(lb_cl_phi_56, 10, 2);
	gbge_cl_nsa->addWidget(le_cl_phi_56, 10, 3);

	/*******************************************************/
	/*********************Aperture radius*******************/
	/*******************************************************/
	lb_cl_ar_min = new QLabel(tr("<big>r</big><sub style=""font-size:large"">min</sub> [mrad]"));

	le_cl_ar_min = new QLineEdit;
	le_cl_ar_min->setStatusTip(tr("Minimum aperture radius"));
	le_cl_ar_min->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_ar_min->setText("0.000");

	lb_cl_ar_max = new QLabel(tr("<big>r</big><sub style=""font-size:large"">max</sub> [mrad]"));

	le_cl_ar_max = new QLineEdit;
	le_cl_ar_max->setStatusTip(tr("Maximum aperture radius"));
	le_cl_ar_max->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_cl_ar_max->setText("21.000");

	gbh_cl_ar = new MGroupBox_H(tr("Aperture radius"));
	gbh_cl_ar->addWidget(lb_cl_ar_min, le_cl_ar_min, lb_cl_ar_max, le_cl_ar_max);

	/*******************************************************/
	/*********************Condenser lens********************/
	/*******************************************************/
	auto layout = new QVBoxLayout;
	layout->setSpacing(0);
	layout->setContentsMargins(2, 1, 1, 0);
	layout->addLayout(lyg_cl_sa);
	layout->addWidget(gbge_cl_nsa);
	layout->addWidget(gbh_cl_ar);
	layout->addStretch();

	lens = new QWidget;
	lens->setLayout(layout);
}

void MainWindow::create_objective_lens_tab(QWidget *&lens)
{
	/*******************************************************/
	/*****************Symmetrical aberrations***************/
	/*******************************************************/
	/************3rd order spherical aberration************/
	lb_ol_sa_c_30 = new QLabel(tr(c_nm_to_qba(30, " [mm]")));

	le_ol_sa_c_30 = new QLineEdit;
	le_ol_sa_c_30->setStatusTip(tr("3rd order spherical aberration"));
	le_ol_sa_c_30->setValidator(new QDoubleValidator(-double_lim, double_lim, 5));
	le_ol_sa_c_30->setText("0.001");

	/************5th order spherical aberration************/
	lb_ol_sa_c_50 = new QLabel(tr(c_nm_to_qba(50, " [mm]")));

	le_ol_sa_c_50 = new QLineEdit;
	le_ol_sa_c_50->setStatusTip(tr("5th order spherical aberration"));
	le_ol_sa_c_50->setValidator(new QDoubleValidator(-double_lim, double_lim, 5));
	le_ol_sa_c_50->setText("0.000");

	/*******************Vortex momentum*********************/
	lb_ol_sa_m = new QLabel(tr("m"));
	lb_ol_sa_m->setAlignment(Qt::AlignCenter);

	le_ol_sa_m = new QLineEdit;
	le_ol_sa_m->setStatusTip(tr("Vortex momentum"));
	le_ol_sa_m->setValidator(new QIntValidator(-int_lim, int_lim));
	le_ol_sa_m->setText("0");

	/*************************Defocus***********************/
	lb_ol_sa_c_10 = new QLabel(tr(c_nm_to_qba(10, " [Å]")));

	le_ol_sa_c_10 = new QLineEdit;
	le_ol_sa_c_10->setStatusTip(tr("Defocus [Å]"));
	le_ol_sa_c_10->setValidator(new QDoubleValidator(-double_lim, double_lim, 3));
	le_ol_sa_c_10->setText("0.000");

	pb_ol_sa_c_10_opt = new QPushButton(QIcon(":/images/opt.png"), "");
	pb_ol_sa_c_10_opt->setStatusTip(tr("Scherzer defocus"));
	pb_ol_sa_c_10_opt->setMaximumSize(22, 22);

	lb_ol_sa_c_10_zero = new QLabel(tr("Zero focus"));

	cb_ol_sa_c_10_zero = new QComboBox;
	cb_ol_sa_c_10_zero->setStatusTip(tr("Zero reference focus"));
	cb_ol_sa_c_10_zero->addItem(tr("First atom"), QVariant(mt::eZDT_First));
	cb_ol_sa_c_10_zero->setItemData(0, tr("Zero focus located at the first atom"), Qt::ToolTipRole);
	cb_ol_sa_c_10_zero->addItem(tr("Half thick."), QVariant(mt::eZDT_Middle));
	cb_ol_sa_c_10_zero->setItemData(1, tr("Zero focus located at half thickness"), Qt::ToolTipRole);
	cb_ol_sa_c_10_zero->addItem(tr("Last atom"), QVariant(mt::eZDT_Last));
	cb_ol_sa_c_10_zero->setItemData(2, tr("Zero focus located at the last atom"), Qt::ToolTipRole);
	cb_ol_sa_c_10_zero->addItem(tr("User def."), QVariant(mt::eZDT_User_Define));
	cb_ol_sa_c_10_zero->setItemData(3, tr("Zero focus located at user define plane"), Qt::ToolTipRole);
	lb_ol_sa_c_10_z = new QLabel(tr("Plane [Å]"));

	le_ol_sa_c_10_z = new QLineEdit;
	le_ol_sa_c_10_z->setStatusTip(tr("Set focus plane"));
	le_ol_sa_c_10_z->setValidator(new QDoubleValidator(-double_lim, double_lim, 3));
	le_ol_sa_c_10_z->setText("0.000");

	/*******************************************************/
	lb_ol_sa_delta = new QLabel(tr("<span style=""font-family:symbol;font-size:large"">D</span> [Å]"));
	lb_ol_sa_delta->setAlignment(Qt::AlignCenter);

	le_ol_sa_delta = new QLineEdit;
	le_ol_sa_delta->setStatusTip(tr("Defocus spread"));
	le_ol_sa_delta->setValidator(new QDoubleValidator(0, std::numeric_limits<double>::max(), 3));
	le_ol_sa_delta->setText("20.0");

	/*******************************************************/
	lb_ol_sa_ndelta = new QLabel(tr("N<sub style=""font-family:symbol;font-size:large"">D</sub>"));
	lb_ol_sa_ndelta->setAlignment(Qt::AlignCenter);

	cb_ol_sa_ndelta = new QComboBox;
	cb_ol_sa_ndelta->setStatusTip(tr("Number of integration points"));
	for(auto ik=2; ik<=128; ik++)
	{
		cb_ol_sa_ndelta->addItem(QString::number(ik), QVariant(ik));
	}
	cb_ol_sa_ndelta->setCurrentIndex(6);

	/*******************************************************/
	auto lyg_ol_sa = new QGridLayout;
	lyg_ol_sa->setHorizontalSpacing(3);
	lyg_ol_sa->setVerticalSpacing(0);
	lyg_ol_sa->setContentsMargins(2, 1, 1, 0);

	lyg_ol_sa->addWidget(lb_ol_sa_c_30, 0, 0);
	lyg_ol_sa->addWidget(le_ol_sa_c_30, 0, 1);
	lyg_ol_sa->addWidget(lb_ol_sa_c_50, 0, 2);
	lyg_ol_sa->addWidget(le_ol_sa_c_50, 0, 3);

	lyg_ol_sa->addWidget(lb_ol_sa_m, 1, 0);
	lyg_ol_sa->addWidget(le_ol_sa_m, 1, 1);
	lyg_ol_sa->addWidget(lb_ol_sa_c_10, 1, 2);

	auto lyh_c_10_opt = new QHBoxLayout;
	lyh_c_10_opt->setSpacing(0);
	lyh_c_10_opt->addWidget(le_ol_sa_c_10);
	lyh_c_10_opt->addWidget(pb_ol_sa_c_10_opt);

	lyg_ol_sa->addLayout(lyh_c_10_opt, 1, 3);

	lyg_ol_sa->addWidget(lb_ol_sa_c_10_zero, 2, 0);
	lyg_ol_sa->addWidget(cb_ol_sa_c_10_zero, 2, 1);
	lyg_ol_sa->addWidget(lb_ol_sa_c_10_z, 2, 2);
	lyg_ol_sa->addWidget(le_ol_sa_c_10_z, 2, 3);

	lyg_ol_sa->addWidget(lb_ol_sa_delta, 3, 0);
	lyg_ol_sa->addWidget(le_ol_sa_delta, 3, 1);
	lyg_ol_sa->addWidget(lb_ol_sa_ndelta, 3, 2);
	lyg_ol_sa->addWidget(cb_ol_sa_ndelta, 3, 3);

	connect(pb_ol_sa_c_10_opt, SIGNAL(released()), SLOT(pb_ol_sa_c_10_opt_released()));
	pb_ol_sa_c_10_opt->click();

	connect(cb_ol_sa_c_10_zero, SIGNAL(currentIndexChanged(int)), SLOT(cb_ol_sa_c_10_zero_currentIndexChanged(int)));
	auto idx = cb_ol_sa_c_10_zero->findData(mt::eZDT_Last);
	cb_ol_sa_c_10_zero->setCurrentIndex(idx);
	cb_ol_sa_c_10_zero_currentIndexChanged(idx);

	/*******************************************************/
	/***************Nonsymmetrical aberrations**************/
	/*******************************************************/

	/*******************2-fold astigmatism******************/
	lb_ol_c_12 = new QLabel(tr(c_nm_to_qba(12, " [Å]")));

	le_ol_c_12 = new QLineEdit;
	le_ol_c_12->setStatusTip(tr("2-fold astigmatism"));
	le_ol_c_12->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_12->setText("0.000");

	lb_ol_phi_12 = new QLabel(tr(phi_nm_to_qba(12, " [°]")));

	le_ol_phi_12 = new QLineEdit;
	le_ol_phi_12->setStatusTip(tr("Azimuthal angle of 2-fold astigmatism"));
	le_ol_phi_12->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_12->setText("0.000");

	/*******************************************************/
	/***************************c2**************************/

	/***********************Axial coma**********************/
	lb_ol_c_21 = new QLabel(tr(c_nm_to_qba(21, " [Å]")));

	le_ol_c_21 = new QLineEdit;
	le_ol_c_21->setStatusTip(tr("Axial coma"));
	le_ol_c_21->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_21->setText("0.000");

	lb_ol_phi_21 = new QLabel(tr(phi_nm_to_qba(21, " [°]")));

	le_ol_phi_21 = new QLineEdit;
	le_ol_phi_21->setStatusTip(tr("Azimuthal angle of axial coma"));
	le_ol_phi_21->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_21->setText("0.000");

	/*******************3-fold astigmatism******************/
	lb_ol_c_23 = new QLabel(tr(c_nm_to_qba(23, " [Å]")));

	le_ol_c_23 = new QLineEdit;
	le_ol_c_23->setStatusTip(tr("3-fold astigmatism"));
	le_ol_c_23->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_23->setText("0.000");

	lb_ol_phi_23 = new QLabel(tr(phi_nm_to_qba(23, " [°]")));

	le_ol_phi_23 = new QLineEdit;
	le_ol_phi_23->setStatusTip(tr("Azimuthal angle of 3-fold astigmatism"));
	le_ol_phi_23->setValidator(new QDoubleValidator(360, 360, 3));
	le_ol_phi_23->setText("0.000");

	/*******************************************************/
	/***************************c3**************************/

	/****************Axial star aberration******************/
	lb_ol_c_32 = new QLabel(tr(c_nm_to_qba(32, " [Å]")));

	le_ol_c_32 = new QLineEdit;
	le_ol_c_32->setStatusTip(tr("Axial star aberration"));
	le_ol_c_32->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_32->setText("0.000");

	lb_ol_phi_32 = new QLabel(tr(phi_nm_to_qba(32, " [°]")));

	le_ol_phi_32 = new QLineEdit;
	le_ol_phi_32->setStatusTip(tr("Azimuthal angle of axial star aberration"));
	le_ol_phi_32->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_32->setText("0.000");

	/*******************4-fold astigmatism******************/
	lb_ol_c_34 = new QLabel(tr(c_nm_to_qba(34, " [Å]")));

	le_ol_c_34 = new QLineEdit;
	le_ol_c_34->setStatusTip(tr("4-fold astigmatism"));
	le_ol_c_34->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_34->setText("0.000");

	lb_ol_phi_34 = new QLabel(tr(phi_nm_to_qba(34, " [°]")));

	le_ol_phi_34 = new QLineEdit;
	le_ol_phi_34->setStatusTip(tr("Azimuthal angle of 4-fold astigmatism"));
	le_ol_phi_34->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_34->setText("0.000");

	/*******************************************************/
	/***************************c4**************************/

	/*****************4th order axial coma******************/
	lb_ol_c_41 = new QLabel(tr(c_nm_to_qba(41, " [Å]")));

	le_ol_c_41 = new QLineEdit;
	le_ol_c_41->setStatusTip(tr("4th order axial coma"));
	le_ol_c_41->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_41->setText("0.000");

	lb_ol_phi_41 = new QLabel(tr(phi_nm_to_qba(41, " [°]")));

	le_ol_phi_41 = new QLineEdit;
	le_ol_phi_41->setStatusTip(tr("Azimuthal angle of 4th order axial coma"));
	le_ol_phi_41->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_41->setText("0.000");

	/*******************3-lobe aberration******************/
	lb_ol_c_43 = new QLabel(tr(c_nm_to_qba(43, " [Å]")));

	le_ol_c_43 = new QLineEdit;
	le_ol_c_43->setStatusTip(tr("3-lobe aberration"));
	le_ol_c_43->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_43->setText("0.000");

	lb_ol_phi_43 = new QLabel(tr(phi_nm_to_qba(43, " [°]")));

	le_ol_phi_43 = new QLineEdit;
	le_ol_phi_43->setStatusTip(tr("Azimuthal angle of 3-lobe aberration"));
	le_ol_phi_43->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_43->setText("0.000");

	/*******************5-fold astigmatism******************/
	lb_ol_c_45 = new QLabel(tr(c_nm_to_qba(45, " [Å]")));

	le_ol_c_45 = new QLineEdit;
	le_ol_c_45->setStatusTip(tr("5-fold astigmatism"));
	le_ol_c_45->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_45->setText("0.000");

	lb_ol_phi_45 = new QLabel(tr(phi_nm_to_qba(45, " [°]")));

	le_ol_phi_45 = new QLineEdit;
	le_ol_phi_45->setStatusTip(tr("Azimuthal angle of 5-fold astigmatism"));
	le_ol_phi_45->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_45->setText("0.000");

	/*******************************************************/
	/***************************c5**************************/

	/*************5th order axial star aberration***********/
	lb_ol_c_52 = new QLabel(tr(c_nm_to_qba(52, " [Å]")));

	le_ol_c_52 = new QLineEdit;
	le_ol_c_52->setStatusTip(tr("5th order axial star aberration"));
	le_ol_c_52->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_52->setText("0.000");

	lb_ol_phi_52 = new QLabel(tr(phi_nm_to_qba(52, " [°]")));

	le_ol_phi_52 = new QLineEdit;
	le_ol_phi_52->setStatusTip(tr("Azimuthal angle of 5th order axial star aberration"));
	le_ol_phi_52->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_52->setText("0.000");

	/**************5th order rosette aberration*************/
	lb_ol_c_54 = new QLabel(tr(c_nm_to_qba(54, " [Å]")));

	le_ol_c_54 = new QLineEdit;
	le_ol_c_54->setStatusTip(tr("5th order rosette aberration"));
	le_ol_c_54->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_54->setText("0.000");

	lb_ol_phi_54 = new QLabel(tr(phi_nm_to_qba(54, " [°]")));

	le_ol_phi_54 = new QLineEdit;
	le_ol_phi_54->setStatusTip(tr("Azimuthal angle of 5th order rosette aberration"));
	le_ol_phi_54->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_54->setText("0.000");

	/*******************6-fold astigmatism******************/
	lb_ol_c_56 = new QLabel(tr(c_nm_to_qba(56, " [Å]")));

	le_ol_c_56 = new QLineEdit;
	le_ol_c_56->setStatusTip(tr("6-fold astigmatism"));
	le_ol_c_56->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_c_56->setText("0.000");

	lb_ol_phi_56 = new QLabel(tr(phi_nm_to_qba(56, " [°]")));

	le_ol_phi_56 = new QLineEdit;
	le_ol_phi_56->setStatusTip(tr("Azimuthal angle of 6-fold astigmatism"));
	le_ol_phi_56->setValidator(new QDoubleValidator(-360, 360, 3));
	le_ol_phi_56->setText("0.000");

	gbge_ol_nsa = new MGroupBox_GE(tr("Nonsymmetrical aberrations"));

	gbge_ol_nsa->addWidget(lb_ol_c_12, 0, 0);
	gbge_ol_nsa->addWidget(le_ol_c_12, 0, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_12, 0, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_12, 0, 3);

	gbge_ol_nsa->addWidget(lb_ol_c_21, 1, 0);
	gbge_ol_nsa->addWidget(le_ol_c_21, 1, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_21, 1, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_21, 1, 3);
	gbge_ol_nsa->addWidget(lb_ol_c_23, 2, 0);
	gbge_ol_nsa->addWidget(le_ol_c_23, 2, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_23, 2, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_23, 2, 3);

	gbge_ol_nsa->addWidget(lb_ol_c_32, 3, 0);
	gbge_ol_nsa->addWidget(le_ol_c_32, 3, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_32, 3, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_32, 3, 3);
	gbge_ol_nsa->addWidget(lb_ol_c_34, 4, 0);
	gbge_ol_nsa->addWidget(le_ol_c_34, 4, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_34, 4, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_34, 4, 3);

	gbge_ol_nsa->addWidget(lb_ol_c_41, 5, 0);
	gbge_ol_nsa->addWidget(le_ol_c_41, 5, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_41, 5, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_41, 5, 3);
	gbge_ol_nsa->addWidget(lb_ol_c_43, 6, 0);
	gbge_ol_nsa->addWidget(le_ol_c_43, 6, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_43, 6, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_43, 6, 3);
	gbge_ol_nsa->addWidget(lb_ol_c_45, 7, 0);
	gbge_ol_nsa->addWidget(le_ol_c_45, 7, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_45, 7, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_45, 7, 3);

	gbge_ol_nsa->addWidget(lb_ol_c_52, 8, 0);
	gbge_ol_nsa->addWidget(le_ol_c_52, 8, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_52, 8, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_52, 8, 3);
	gbge_ol_nsa->addWidget(lb_ol_c_54, 9, 0);
	gbge_ol_nsa->addWidget(le_ol_c_54, 9, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_54, 9, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_54, 9, 3);
	gbge_ol_nsa->addWidget(lb_ol_c_56, 10, 0);
	gbge_ol_nsa->addWidget(le_ol_c_56, 10, 1);
	gbge_ol_nsa->addWidget(lb_ol_phi_56, 10, 2);
	gbge_ol_nsa->addWidget(le_ol_phi_56, 10, 3);

	/*******************************************************/
	/*********************Aperture radius*******************/
	/*******************************************************/
	lb_ol_ar_min = new QLabel(tr("<big>r</big><sub style=""font-size:large"">min</sub> [mrad]"));

	le_ol_ar_min = new QLineEdit;
	le_ol_ar_min->setStatusTip(tr("Minimum aperture radius"));
	le_ol_ar_min->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_ar_min->setText("0.000");

	lb_ol_ar_max = new QLabel(tr("<big>r</big><sub style=""font-size:large"">max</sub> [mrad]"));

	le_ol_ar_max = new QLineEdit;
	le_ol_ar_max->setStatusTip(tr("Maximum aperture radius"));
	le_ol_ar_max->setValidator(new QDoubleValidator(0, double_lim, 3));
	le_ol_ar_max->setText("0.000");

	gbh_ol_ar = new MGroupBox_H(tr("Aperture radius"));
	gbh_ol_ar->addWidget(lb_ol_ar_min, le_ol_ar_min, lb_ol_ar_max, le_ol_ar_max);

	/*******************************************************/
	/*********************Objective lens********************/
	/*******************************************************/
	auto layout = new QVBoxLayout;
	layout->setSpacing(0);
	layout->setContentsMargins(2, 1, 1, 0);
	layout->addLayout(lyg_ol_sa);
	layout->addWidget(gbge_ol_nsa);
	layout->addWidget(gbh_ol_ar);
	layout->addStretch();

	lens = new QWidget;
	lens->setLayout(layout);
}

/*******************************************************/
void MainWindow::create_general_dock_widget()
{
	auto fr_general = new MFrame_V;
	fr_general->layout->setSpacing(1);
	fr_general->setStyleSheet("QFrame { background: "+ QColor("#D8D8D8").name() +" }");

	/*******************************************************/
	/************************device*************************/
	/*******************************************************/
	cb_device = new QComboBox;
	cb_device->setStatusTip(tr("Device calculation"));

	cb_device->setMinimumWidth(80);
	cb_device->addItem(tr("CPU"), mt::e_host);
	cb_device->setItemData(0, tr("CPU calculation"), Qt::ToolTipRole);
	cb_device->addItem(tr("GPU"), mt::e_device);
	cb_device->setItemData(1, tr("GPU calculation"), Qt::ToolTipRole);

	/*******************************************************/
	/****************floating point calculation*************/
	/*******************************************************/
	lb_precision = new QLabel(tr("Precision"));

	cb_precision = new QComboBox;
	cb_precision->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
	cb_precision->setStatusTip(tr("Floating-point precision"));
	cb_precision->addItem(tr("Single"), enum_to_QVar(mt::eP_float));
	cb_precision->setItemData(0, tr("Single precision (32 bits)"), Qt::ToolTipRole);
	cb_precision->addItem(tr("Double"), enum_to_QVar(mt::eP_double));
	cb_precision->setItemData(1, tr("Double precision (64 bits)"), Qt::ToolTipRole);
	cb_precision->setCurrentIndex(0);

	/*******************************************************/
	/*****************Number of threads*********************/
	/*******************************************************/
	lb_nthreads = new QLabel(tr(x_sub_to_qba("N", "threads", "x-large")));
	lb_nthreads->setAlignment(Qt::AlignCenter);

	sb_nthreads = new QSpinBox;
	sb_nthreads->setStatusTip(tr("Number of threads"));
	auto nthreads = QThread::idealThreadCount();
	sb_nthreads->setRange(1, nthreads);
	sb_nthreads->setSingleStep(1);
	sb_nthreads->setValue(nthreads);

	/*******************************************************/
	/***********************GPU card************************/
	/*******************************************************/
	lb_gpu_card = new QLabel(tr("Graphic cards"));
	lb_gpu_card->setAlignment(Qt::AlignCenter);

	cb_gpu_card = new QComboBox;
	cb_gpu_card->setStatusTip(tr("Available graphic cards"));
	cb_gpu_card->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Minimum);

	  std::vector<mt::Device_Properties> device_properties;
	  mt::get_device_properties(device_properties);
	  for(auto id=0; id<device_properties.size(); id++)
	  {
		cb_gpu_card->addItem(QString::fromStdString(device_properties[id].name), QVariant(device_properties[id].id));
	  }

	/*******************************************************/
	gbg_device = new MGroupBox_G(tr("Device"));
	gbg_device->layout->setVerticalSpacing(2);
	gbg_device->layout->setContentsMargins(2, 0, 1, 0);

	gbg_device->addWidget(cb_device, 0, 0);

	auto lyh_precision= new QHBoxLayout;
	lyh_precision->setSpacing(2);
	lyh_precision->setContentsMargins(0, 0, 0, 1);
	lyh_precision->addWidget(lb_precision);
	lyh_precision->addWidget(cb_precision);

	gbg_device->layout->addLayout(lyh_precision, 0, 1);

	gbg_device->addWidget(lb_nthreads, 1, 0);
	gbg_device->addWidget(sb_nthreads, 1, 1);
	gbg_device->addWidget(lb_gpu_card, 1, 0);
	gbg_device->addWidget(cb_gpu_card, 1, 1);

	fr_general->layout->addWidget(gbg_device);

	/*******************************************************/
	connect(cb_device, SIGNAL(currentIndexChanged(int)), SLOT(cb_device_currentIndexChanged(int)));

	auto idx = (mt::is_gpu_available())?1:0;
	idx = 0;
	cb_device->setCurrentIndex(idx);
	cb_device_currentIndexChanged(idx);

	/*******************************************************/
	/*********Electron-Specimen interation model************/
	/*******************************************************/
	cb_elec_spec_int_model = new QComboBox;
	cb_elec_spec_int_model->setStatusTip(tr("Electron-Specimen interaction model"));
	cb_elec_spec_int_model->addItem(tr("Multislice approx."), mt::eESIM_Multislice);
	cb_elec_spec_int_model->setItemData(0, tr("Multislice approximation"), Qt::ToolTipRole);
	cb_elec_spec_int_model->addItem(tr("Phase object approx."), mt::eESIM_Phase_Object);
	cb_elec_spec_int_model->setItemData(1, tr("Phase object approximation"), Qt::ToolTipRole);
	cb_elec_spec_int_model->addItem(tr("Weak phase object approx."), mt::eESIM_Weak_Phase_Object);
	cb_elec_spec_int_model->setItemData(2, tr("Weak phase object approximation"), Qt::ToolTipRole);

	/*******************************************************/
	/***********Parameterization of the fe(g)***************/
	/*******************************************************/
	lb_potential_type = new QLabel(tr("fₑ(g)"));
	lb_potential_type->setAlignment(Qt::AlignCenter);

	cb_potential_type = new QComboBox;
	cb_potential_type->setStatusTip(tr("Parameterization of the electron scattering factors"));
	cb_potential_type->addItem(tr("Doyle et al. 0 < g < 4")+rAgs, enum_to_QVar(mt::ePT_Doyle_0_4));
	cb_potential_type->setItemData(0, tr("fₑ(g) was fitted for scattering angles 0 < g < 4")+rAgs, Qt::ToolTipRole);
	cb_potential_type->addItem(tr("Peng et al. 0 < g < 4")+rAgs, enum_to_QVar(mt::ePT_Peng_0_4));
	cb_potential_type->setItemData(1, tr("fₑ(g) was fitted for scattering angles 0 < g < 4")+rAgs, Qt::ToolTipRole);
	cb_potential_type->addItem(tr("Peng et al. 0 < g < 12")+rAgs, enum_to_QVar(mt::ePT_Peng_0_12));
	cb_potential_type->setItemData(2, tr("fₑ(g) was fitted for scattering angles 0 < g < 12")+rAgs, Qt::ToolTipRole);
	cb_potential_type->addItem(tr("Kirkland 0 < g < 12")+rAgs, enum_to_QVar(mt::ePT_Kirkland_0_12));
	cb_potential_type->setItemData(3, tr("fₑ(g) was fitted for scattering angles 0 < g < 12")+rAgs, Qt::ToolTipRole);
	cb_potential_type->addItem(tr("Weickenmeier et al. 0 < g < 12")+rAgs, enum_to_QVar(mt::ePT_Weickenmeier_0_12));
	cb_potential_type->setItemData(4, tr("fₑ(g) was fitted for scattering angles 0 < g < 12")+rAgs, Qt::ToolTipRole);
	cb_potential_type->addItem(tr("Lobato et al. 0 < g < 12")+rAgs, enum_to_QVar(mt::ePT_Lobato_0_12));
	cb_potential_type->setItemData(5, tr("fₑ(g) was fitted for scattering angles 0 < g < 12")+rAgs, Qt::ToolTipRole);
	cb_potential_type->setCurrentIndex(5);

	/*******************************************************/
	gbg_elec_spec_int_model = new MGroupBox_G(elec + tr(" - Specimen interaction"));
	gbg_elec_spec_int_model->layout->setVerticalSpacing(1);
	gbg_elec_spec_int_model->layout->setHorizontalSpacing(4);

	gbg_elec_spec_int_model->addWidget(cb_elec_spec_int_model, 0, 0, 1, 2);
	gbg_elec_spec_int_model->addWidget(lb_potential_type, 1, 0);
	gbg_elec_spec_int_model->addWidget(cb_potential_type, 1, 1);

	fr_general->layout->addWidget(gbg_elec_spec_int_model);

	/*******************************************************/
	/************Electron-Phonon interation model***********/
	/*******************************************************/
	cb_pn_model = new QComboBox;
	cb_pn_model->setStatusTip(tr("Electron-Phonon interaction model"));
	cb_pn_model->addItem(tr("Still atom approx."), enum_to_QVar(mt::ePM_Still_Atom));
	cb_pn_model->setItemData(0, tr("Still atom approximation"), Qt::ToolTipRole);
	cb_pn_model->addItem(tr("Absorptive pot. approx."), enum_to_QVar(mt::ePM_Absorptive_Model));
	cb_pn_model->setItemData(1, tr("Absorptive potential approximation"), Qt::ToolTipRole);
	cb_pn_model->addItem(tr("Frozen phonon approx."), enum_to_QVar(mt::ePM_Frozen_Phonon));
	cb_pn_model->setItemData(2, tr("Frozen phonon approximation"), Qt::ToolTipRole);
	cb_pn_model->setCurrentIndex(2);

	disable_item_QComboBox(cb_pn_model, 1);

	/*******************************************************/
	ckb_pn_coh_contrib = new QCheckBox(tr("Coh. contribution"));
	ckb_pn_coh_contrib->setStatusTip(tr("Calculate coherent contribution"));

	/*******************************************************/
	ckb_pn_single_conf = new QCheckBox(tr("Sng. configuration"));
	ckb_pn_single_conf->setStatusTip(tr("Frozen phonon for a specific configuration"));

	/*******************************************************/
	lb_pn_nconf = new QLabel(tr(x_sub_to_qba("N", "conf.", "x-large")));

	le_pn_nconf = new QLineEdit;
	le_pn_nconf->setStatusTip(tr("Number of configurations"));
	le_pn_nconf->setValidator(new QIntValidator(0, std::numeric_limits<int>::max()));
	le_pn_nconf->setText("5");

	/*******************************************************/
	lb_pn_dim = new QLabel(tr("Vibration"));

	cb_pn_dim = new QComboBox;
	cb_pn_dim->setStatusTip(tr("Direction of the atomic vibrations"));
	cb_pn_dim->addItem("x", QVariant(100));
	cb_pn_dim->setItemData(0, "x vibration", Qt::ToolTipRole);
	cb_pn_dim->addItem("y", QVariant(10));
	cb_pn_dim->setItemData(1, "y vibration", Qt::ToolTipRole);
	cb_pn_dim->addItem("z", QVariant(1));
	cb_pn_dim->setItemData(2, "z vibration", Qt::ToolTipRole);
	cb_pn_dim->addItem("xy", QVariant(110));
	cb_pn_dim->setItemData(3, "x-y vibration", Qt::ToolTipRole);
	cb_pn_dim->addItem("xz", QVariant(101));
	cb_pn_dim->setItemData(4, "x-z vibration", Qt::ToolTipRole);
	cb_pn_dim->addItem("yz", QVariant(11));
	cb_pn_dim->setItemData(5, "x-z vibration", Qt::ToolTipRole);
	cb_pn_dim->addItem("xyz", QVariant(111));
	cb_pn_dim->setItemData(6, "x-y-z vibration", Qt::ToolTipRole);
	cb_pn_dim->setCurrentIndex(cb_pn_dim->findData(110));

	/*******************************************************/
	lb_pn_seed = new QLabel(tr("Rnd seed"));

	le_pn_seed = new QLineEdit;
	le_pn_seed->setStatusTip(tr("Random generator seed"));
	le_pn_seed->setValidator(new QIntValidator(0, int_lim));
	le_pn_seed->setText("300183");

	/*******************************************************/
	ckb_pn_auto_seed = new QCheckBox(tr("auto"));
	ckb_pn_auto_seed->setStatusTip(tr("Automatic cgpu_rand seed"));

	/*******************************************************/
	auto lyg_frozen_phonon = new QGridLayout;
	lyg_frozen_phonon->setVerticalSpacing(0);
	lyg_frozen_phonon->setHorizontalSpacing(2);
	lyg_frozen_phonon->setContentsMargins(0, 0, 0, 0);

	lyg_frozen_phonon->addWidget(ckb_pn_coh_contrib, 0, 0, 1, 2);
	lyg_frozen_phonon->addWidget(ckb_pn_single_conf, 0, 2, 1, 2);

	lyg_frozen_phonon->addWidget(lb_pn_nconf, 1, 0);
	lyg_frozen_phonon->addWidget(le_pn_nconf, 1, 1);
	lyg_frozen_phonon->addWidget(lb_pn_dim, 1, 2);
	lyg_frozen_phonon->addWidget(cb_pn_dim, 1, 3);

	lyg_frozen_phonon->addWidget(lb_pn_seed, 2, 0);

	lyg_frozen_phonon->addWidget(le_pn_seed, 2, 1);
	lyg_frozen_phonon->addWidget(ckb_pn_auto_seed, 2, 2);

	wg_electron_phonon = new QWidget;
	wg_electron_phonon->setLayout(lyg_frozen_phonon);

	/*******************************************************/
	gbg_elec_phonon_int_model = new MGroupBox_G(elec + tr(" - Phonon interaction"));
	gbg_elec_phonon_int_model->layout->setVerticalSpacing(1);

	gbg_elec_phonon_int_model->addWidget(cb_pn_model, 0, 0);
	gbg_elec_phonon_int_model->addWidget(wg_electron_phonon, 1, 0);

	fr_general->layout->addWidget(gbg_elec_phonon_int_model);
	fr_general->layout->addStretch(1);

	/*******************************************************/
	connect(cb_pn_model, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
			[this](int index)
	{
		auto pn_model = QVar_to_enum<mt::ePhonon_Model>(cb_pn_model->itemData(index));
		auto bb_frozen_phonon = pn_model == mt::ePM_Frozen_Phonon;
		wg_electron_phonon->setEnabled(bb_frozen_phonon);
	});

	connect(ckb_pn_single_conf, &QCheckBox::stateChanged,
			[this](int state)
	{
		if(state==Qt::Checked)
		{
			lb_pn_nconf->setText(tr("K<sub style=""font-size:x-large"">conf.</sub>"));
			le_pn_nconf->setStatusTip(tr("Specific configuration"));
		}
		else
		{
			lb_pn_nconf->setText(tr(x_sub_to_qba("N", "conf.", "x-large")));
			le_pn_nconf->setStatusTip(tr("Number of configurations"));
		}
	});

	connect(ckb_pn_auto_seed, &QCheckBox::stateChanged,
			[this](int state)
	{
		auto bb = state==Qt::Checked;
		le_pn_seed->setEnabled(!bb);
		if(bb)
		{
			le_pn_seed->setText(QString::number(qrand()));
		}

	});
	/*******************************************************/
	dw_general = new QDockWidget(tr("General settings"), this);
	dw_general->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	dw_general->setFeatures(QDockWidget::NoDockWidgetFeatures);
	dw_general->setWidget(fr_general);
	addDockWidget(Qt::LeftDockWidgetArea, dw_general);
}

void MainWindow::create_specimen_dock_widget()
{
	//////////////////////////////////////////////
	auto fr_specimen = new MFrame_V;
	fr_specimen->layout->setSpacing(1);
	fr_specimen->setStyleSheet("QFrame { background: "+ QColor("#D8D8D8").name() +" }");

	lb_spec_file = new QLabel(tr("Filename"));

	le_spec_file = new QLineEdit;
	le_spec_file->setStatusTip(tr("Specimen filename"));
	le_spec_file->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);

	set_non_editable_QLineEdit(le_spec_file);

	pb_spec_load = new QPushButton(QIcon(":/images/open.png"), "");
	pb_spec_load->setStatusTip(tr("Load specimen"));
	pb_spec_load->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	pb_spec_show = new QPushButton(QIcon(":/images/view.png"), "");
	pb_spec_show->setStatusTip(tr("Show specimen"));
	pb_spec_show->setMaximumSize(22, 22);
	pb_spec_show->setEnabled(false);

	auto lyh_spec_file = new QHBoxLayout;
	lyh_spec_file->setSpacing(1);
	lyh_spec_file->setContentsMargins(2, 0, 1, 0);
	lyh_spec_file->addWidget(lb_spec_file);
	lyh_spec_file->addSpacing(2);
	lyh_spec_file->addWidget(le_spec_file);
	lyh_spec_file->addWidget(pb_spec_load);
	lyh_spec_file->addWidget(pb_spec_show);

	fr_specimen->layout->addLayout(lyh_spec_file);

	/*******************************************************/
	/**********************Specimen info********************/
	/*******************************************************/
	lb_spec_n_atom_types = new QLabel(tr(n_atoms_types_qba(0, 0)));
	lb_spec_n_atom_types->setFrameStyle(QFrame::Box | QFrame::Sunken);
	lb_spec_n_atom_types->setAlignment(Qt::AlignCenter);
	lb_spec_n_atom_types->setMinimumHeight(25);

	/*******************************************************/
	lb_spec_lx = new QLabel(tr(x_sub_to_qba("L", "x")));
	lb_spec_lx->setAlignment(Qt::AlignCenter);
	lb_spec_lx->setMinimumWidth(20);

	le_spec_lx = new QLineEdit;
	le_spec_lx->setStatusTip(tr("Simulation box size along x"));
	le_spec_lx->setValidator(new QDoubleValidator(0, double_lim, 4));
	le_spec_lx->setText("0.000");

	lb_spec_ly = new QLabel(tr(x_sub_to_qba("L", "y")));
	lb_spec_ly->setAlignment(Qt::AlignCenter);
	lb_spec_ly->setMinimumWidth(20);

	le_spec_ly = new QLineEdit;
	le_spec_ly->setStatusTip(tr("Simulation box size along y"));
	le_spec_ly->setValidator(new QDoubleValidator(0, double_lim, 4));
	le_spec_ly->setText("0.000");

	pb_spec_recenter = new QPushButton("Recenter");
	pb_spec_recenter->setStatusTip(tr("Recenter specimen in x-y directions"));
	pb_spec_recenter->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	auto lyt_dim = new QHBoxLayout;
	lyt_dim->setSpacing(1);
	lyt_dim->setContentsMargins(0, 0, 0, 0);

	lyt_dim->addWidget(lb_spec_lx);
	lyt_dim->addWidget(le_spec_lx);
	lyt_dim->addWidget(lb_spec_ly);
	lyt_dim->addWidget(le_spec_ly);
	lyt_dim->addWidget(pb_spec_recenter);

	/*******************************************************/
	gbg_spec_info = new MGroupBox_G(tr("Specimen info"));
	gbg_spec_info->layout->setVerticalSpacing(1);
	gbg_spec_info->layout->setHorizontalSpacing(2);
	gbg_spec_info->layout->setContentsMargins(2, 0, 1, 0);

	gbg_spec_info->addWidget(lb_spec_n_atom_types, 0, 0);
	gbg_spec_info->layout->addLayout(lyt_dim, 1, 0);

	fr_specimen->layout->addWidget(gbg_spec_info);

	/*******************************************************/
	/*********************** Rotation **********************/
	/*******************************************************/
	lb_spec_rot_theta = new QLabel(tr("<span style=""font-family:symbol;font-size:large"">q</span> [°]"));
	lb_spec_rot_theta->setAlignment(Qt::AlignCenter);

	le_spec_rot_theta = new QLineEdit;
	le_spec_rot_theta->setStatusTip(tr("Rotation angle"));
	le_spec_rot_theta->setValidator(new QDoubleValidator(-360, 360, 3));
	le_spec_rot_theta->setText("0.000");

	lb_spec_rot_u = new QLabel(tr("û"));
	lb_spec_rot_u->setAlignment(Qt::AlignCenter);
	lb_spec_rot_u->setMinimumWidth(15);

	le_spec_rot_ux = new QLineEdit;
	le_spec_rot_ux->setStatusTip(tr("x component of the rotation vector û"));
	le_spec_rot_ux->setValidator(new QDoubleValidator(-double_lim, double_lim, 1));
	le_spec_rot_ux->setText("0.0");
	le_spec_rot_ux->setMaximumWidth(27);

	le_spec_rot_uy = new QLineEdit;
	le_spec_rot_uy->setStatusTip(tr("y component of the rotation vector û"));
	le_spec_rot_uy->setValidator(new QDoubleValidator(-double_lim, double_lim, 1));
	le_spec_rot_uy->setText("0.0");
	le_spec_rot_uy->setMaximumWidth(27);

	le_spec_rot_uz = new QLineEdit;
	le_spec_rot_uz->setStatusTip(tr("z component of the rotation vector û"));
	le_spec_rot_uz->setValidator(new QDoubleValidator(-double_lim, double_lim, 1));
	le_spec_rot_uz->setText("1.0");
	le_spec_rot_uz->setMaximumWidth(27);

	auto lyt_u = new QHBoxLayout;
	lyt_u->setSpacing(1);
	lyt_u->setContentsMargins(0, 0, 0, 0);
	lyt_u->addWidget(le_spec_rot_ux);
	lyt_u->addWidget(le_spec_rot_uy);
	lyt_u->addWidget(le_spec_rot_uz);

	/*******************************************************/
	lb_spec_rot_center_type = new QLabel(tr("R. point"));
	lb_spec_rot_center_type->setAlignment(Qt::AlignCenter);

	cb_spec_rot_center_type = new QComboBox;
	cb_spec_rot_center_type->setStatusTip(tr("Rotation point"));
	cb_spec_rot_center_type->addItem(tr("Centroid"), mt::eRPT_geometric_center);
	cb_spec_rot_center_type->setItemData(0, tr("Geometric center"), Qt::ToolTipRole);
	cb_spec_rot_center_type->addItem(tr("User def."), mt::eRPT_User_Define);
	cb_spec_rot_center_type->setItemData(1, tr("User define center"), Qt::ToolTipRole);

	lb_spec_rot_center_p = new QLabel(tr(x_sub_to_qba("P", "0")));
	lb_spec_rot_center_p->setAlignment(Qt::AlignCenter);
	lb_spec_rot_center_p->setMinimumWidth(15);

	le_spec_rot_center_px = new QLineEdit;
	le_spec_rot_center_px->setStatusTip(tr("x component of the rotation center p0"));
	le_spec_rot_center_px->setValidator(new QDoubleValidator(-double_lim, double_lim, 1));
	le_spec_rot_center_px->setText("0.0");
	le_spec_rot_center_px->setMaximumWidth(27);

	le_spec_rot_center_py = new QLineEdit;
	le_spec_rot_center_py->setStatusTip(tr("y component of the rotation center p0"));
	le_spec_rot_center_py->setValidator(new QDoubleValidator(-double_lim, double_lim, 1));
	le_spec_rot_center_py->setText("0.0");
	le_spec_rot_center_py->setMaximumWidth(27);

	le_spec_rot_center_pz = new QLineEdit;
	le_spec_rot_center_pz->setStatusTip(tr("z component of the rotation center p0"));
	le_spec_rot_center_pz->setValidator(new QDoubleValidator(-double_lim, double_lim, 1));
	le_spec_rot_center_pz->setText("0.0");
	le_spec_rot_center_pz->setMaximumWidth(27);

	auto lyt_p0 = new QHBoxLayout;
	lyt_p0->setSpacing(1);
	lyt_p0->setContentsMargins(0, 0, 0, 0);
	lyt_p0->addWidget(le_spec_rot_center_px);
	lyt_p0->addWidget(le_spec_rot_center_py);
	lyt_p0->addWidget(le_spec_rot_center_pz);

	/*******************************************************/
	gbg_spec_rot = new MGroupBox_G(tr("Rotation"));
	gbg_spec_rot->layout->setVerticalSpacing(1);
	gbg_spec_rot->layout->setHorizontalSpacing(2);
	gbg_spec_rot->layout->setContentsMargins(2, 0, 1, 0);

	gbg_spec_rot->addWidget(lb_spec_rot_theta, 0, 0);
	gbg_spec_rot->addWidget(le_spec_rot_theta, 0, 1);
	gbg_spec_rot->addWidget(lb_spec_rot_u, 0, 2);
	gbg_spec_rot->layout->addLayout(lyt_u, 0, 3);
	gbg_spec_rot->addWidget(lb_spec_rot_center_type, 1, 0);
	gbg_spec_rot->addWidget(cb_spec_rot_center_type, 1, 1);
	gbg_spec_rot->addWidget(lb_spec_rot_center_p, 1, 2);
	gbg_spec_rot->layout->addLayout(lyt_p0, 1, 3);

	fr_specimen->layout->addWidget(gbg_spec_rot);

	connect(le_spec_lx, SIGNAL(textChanged(QString)), SLOT(le_spec_lx_ly_textChanged(QString)));
	connect(le_spec_ly, SIGNAL(textChanged(QString)), SLOT(le_spec_lx_ly_textChanged(QString)));

	connect(pb_spec_recenter, SIGNAL(released()), SLOT(pb_spec_recenter_released()));

	connect(cb_spec_rot_center_type, SIGNAL(currentIndexChanged(int)), SLOT(cb_spec_rot_center_type_currentIndexChanged(int)));
	cb_spec_rot_center_type_currentIndexChanged(cb_spec_rot_center_type->currentIndex());

	/*******************************************************/
	/*********************** Thickness *********************/
	/*******************************************************/
	lb_thk_type = new QLabel(tr("Type"));
	lb_thk_type->setAlignment(Qt::AlignCenter);
	lb_thk_type->setFrameStyle(QFrame::Box | QFrame::Sunken);

	cb_thk_type = new QComboBox;
	cb_thk_type->setStatusTip(tr("Thickness type"));
	cb_thk_type->addItem(tr("whole specimen"), QVariant(mt::eTT_Whole_Spec));
	cb_thk_type->setItemData(0, tr("Whole specimen"), Qt::ToolTipRole);
	cb_thk_type->addItem(tr("by thickness"), QVariant(mt::eTT_Through_Thick));
	cb_thk_type->setItemData(1, tr("By thickness"), Qt::ToolTipRole);
	cb_thk_type->addItem(tr("by slices"), QVariant(mt::eTT_Through_Slices));
	cb_thk_type->setItemData(2, tr("By slices"), Qt::ToolTipRole);

	disable_item_QComboBox(cb_thk_type, 2);

	/***************************************************/
	lb_thk_0 = new QLabel(tr("Thick. min."));
	lb_thk_0->setAlignment(Qt::AlignCenter);
	lb_thk_0->setFrameStyle(QFrame::Box | QFrame::Sunken);

	le_thk_0 = new QLineEdit;
	le_thk_0->setStatusTip(tr("Minimum thickness"));
	le_thk_0->setValidator(new QDoubleValidator(0, double_lim, 4));
	le_thk_0->setText("0.0000");

	/***************************************************/
	lb_thk_d = new QLabel(tr("Thick. step"));
	lb_thk_d->setAlignment(Qt::AlignCenter);
	lb_thk_d->setFrameStyle(QFrame::Box | QFrame::Sunken);

	le_thk_d = new QLineEdit;
	le_thk_d->setStatusTip(tr("Thickness step"));
	le_thk_d->setValidator(new QDoubleValidator(0, double_lim, 4));
	le_thk_d->setText("0.0000");

	/***************************************************/
	lb_thk_e = new QLabel(tr("Thick. max."));
	lb_thk_e->setAlignment(Qt::AlignCenter);
	lb_thk_e->setFrameStyle(QFrame::Box | QFrame::Sunken);

	le_thk_e = new QLineEdit;
	le_thk_e->setStatusTip(tr("Maximum thickness"));
	le_thk_e->setValidator(new QDoubleValidator(0, double_lim, 4));
	le_thk_e->setText("0.0000");

	/*******************************************************/
	gbg_thk = new MGroupBox_G(tr("Thickness"));
	gbg_thk->layout->setVerticalSpacing(1);
	gbg_thk->layout->setHorizontalSpacing(2);
	gbg_thk->layout->setContentsMargins(2, 0, 1, 0);

	gbg_thk->addWidget(lb_thk_type, 0, 0);
	gbg_thk->addWidget(cb_thk_type, 0, 1, 1, 2);
	gbg_thk->addWidget(lb_thk_0, 1, 0);
	gbg_thk->addWidget(lb_thk_d, 1, 1);
	gbg_thk->addWidget(lb_thk_e, 1, 2);
	gbg_thk->addWidget(le_thk_0, 2, 0);
	gbg_thk->addWidget(le_thk_d, 2, 1);
	gbg_thk->addWidget(le_thk_e, 2, 2);

	fr_specimen->layout->addWidget(gbg_thk);

	connect(cb_thk_type, SIGNAL(currentIndexChanged(int)), SLOT(cb_thk_type_currentIndexChanged(int)));

	/*******************************************************/
	/***************** Potential slicing *******************/
	/*******************************************************/
	cb_pot_slic_type = new QComboBox;
	cb_pot_slic_type->setStatusTip(tr("Potential slicing scheme"));
	cb_pot_slic_type->addItem(tr("Slicing by planes"), mt::ePS_Planes);
	cb_pot_slic_type->setItemData(0, tr("Standard potential slicing approach by automatic identification of planes ")+
								   tr("perpedincular to the beam direction"), Qt::ToolTipRole);
	cb_pot_slic_type->addItem(tr("Slicing by dz"), mt::ePS_dz_Proj);
	cb_pot_slic_type->setItemData(1, tr("Standard potential slicing approach using the given slice thickness"), Qt::ToolTipRole);
	cb_pot_slic_type->addItem(tr("Subslicing by dz"), mt::ePS_dz_Sub);
	cb_pot_slic_type->setItemData(2, tr("3d potential subslicing using the given slice thickness"), Qt::ToolTipRole);
	cb_pot_slic_type->addItem(tr("Auto"), mt::ePS_Auto);
	cb_pot_slic_type->setItemData(3, tr("Automatic potential slicing"), Qt::ToolTipRole);

	disable_item_QComboBox(cb_pot_slic_type, 3);

	lb_pot_slic_thick = new QLabel(tr("dz [Å]"));

	le_pot_slic_thick = new QLineEdit;
	le_pot_slic_thick->setStatusTip(tr("Slice thickness"));
	le_pot_slic_thick->setValidator(new QDoubleValidator(0, double_lim, 4));
	le_pot_slic_thick->setText("2.0000");

	pb_pot_slic_view = new QPushButton(QIcon(":/images/view.png"), "");
	pb_pot_slic_view->setStatusTip(tr("Show specimen slicing"));
	pb_pot_slic_view->setMaximumSize(22, 22);
	pb_pot_slic_view->setEnabled(false);

	gbg_potential_slic = new MGroupBox_G(tr("Potential slicing"));
	gbg_potential_slic->setContentsMargins(1, 2, 1, 0);
	gbg_potential_slic->addWidget(cb_pot_slic_type, 0, 0);
	gbg_potential_slic->addWidget(lb_pot_slic_thick, 0, 1);
	gbg_potential_slic->addWidget(le_pot_slic_thick, 0, 2);
	gbg_potential_slic->addWidget(pb_pot_slic_view, 0, 3);

	fr_specimen->layout->addWidget(gbg_potential_slic);

	connect(cb_pot_slic_type, SIGNAL(currentIndexChanged(int)), SLOT(cb_pot_slic_type_currentIndexChanged(int)));
	cb_pot_slic_type->setCurrentIndex(0);
	cb_pot_slic_type_currentIndexChanged(0);

	connect(le_pot_slic_thick, SIGNAL(editingFinished()), SLOT(le_pot_slic_thick_editingFinished()));

	connect(pb_spec_load, SIGNAL(released()), SLOT(pb_spec_load_released()));

	/*******************************************************/
	/********************** xy sampling ********************/
	/*******************************************************/
	lb_samp_nx = new QLabel(tr(x_sub_to_qba("N", "x")));

	le_samp_nx = new QLineEdit;
	le_samp_nx->setStatusTip(tr("Number of pixels in the x direction"));
	le_samp_nx->setValidator(new QIntValidator(32, int_lim));
	le_samp_nx->setText("1024");

	lb_samp_ny = new QLabel(tr(x_sub_to_qba("N", "y")));

	le_samp_ny = new QLineEdit;
	le_samp_ny->setStatusTip(tr("Number of pixels in the y direction"));
	le_samp_ny->setValidator(new QIntValidator(32, int_lim));
	le_samp_ny->setText("1024");

	ckb_samp_bandwidth = new QCheckBox(tr("bwl"));
	ckb_samp_bandwidth->setStatusTip(tr("Symmetrically bandwidth limit the wave function"));

	gbg_samp = new MGroupBox_G(tr("Box sampling"));
	gbg_samp->layout->setVerticalSpacing(1);
	gbg_samp->layout->setHorizontalSpacing(3);

	gbg_samp->addWidget(lb_samp_nx, 0, 0);
	gbg_samp->addWidget(le_samp_nx, 0, 1);
	gbg_samp->addWidget(lb_samp_ny, 0, 2);
	gbg_samp->addWidget(le_samp_ny, 0, 3);
	gbg_samp->addWidget(ckb_samp_bandwidth, 0, 4);

	fr_specimen->layout->addWidget(gbg_samp);
	fr_specimen->layout->addStretch(1);

	/*******************************************************/
	connect(le_samp_nx, &QLineEdit::editingFinished,
	[this]()
	{
		mt::Prime_Num pn;
		auto nx = le_samp_nx->text().toInt();
		nx = pn(nx, mt::eDST_Closest);
		le_samp_nx->setText(QString::number(nx));
	});

	connect(le_samp_ny, &QLineEdit::editingFinished,
	[this]()
	{
		mt::Prime_Num pn;
		auto ny = le_samp_ny->text().toInt();
		ny = pn(ny, mt::eDST_Closest);
		le_samp_ny->setText(QString::number(ny));
	});

	//////////////////////////////////////////////
	dw_specimen = new QDockWidget(tr("Specimen"), this);
	dw_specimen->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	dw_specimen->setFeatures(QDockWidget::NoDockWidgetFeatures);
	dw_specimen->setWidget(fr_specimen);
	addDockWidget(Qt::LeftDockWidgetArea, dw_specimen);
}

void MainWindow::create_experiment_dock_widget()
{
	//////////////////////////////////////////////
	auto *fr_experiment = new MFrame_G;
	fr_experiment->layout->setVerticalSpacing(0);
	fr_experiment->setStyleSheet("QFrame { background: "+ QColor("#D8D8D8").name() +" }");

	lb_sim_type = new QLabel(tr("Select simulation type"));
	lb_sim_type->setAlignment(Qt::AlignCenter);
	//  lb_sim_type->setFrameStyle(QFrame::Box | QFrame::Sunken);

	cb_sim_type = new QComboBox;
	cb_sim_type->setStatusTip(tr("Select simulation type"));

	cb_sim_type->addItem(tr("STEM"), mt::eTEMST_STEM);
	cb_sim_type->setItemData(0, tr("Scanning transmission electron microscopy"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("ISTEM"), mt::eTEMST_ISTEM);
	cb_sim_type->setItemData(1, tr("Imaging scanning transmission electron microscopy"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("CBED"), mt::eTEMST_CBED);
	cb_sim_type->setItemData(2, tr("Convergent beam electron diffraction"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("CBEI"), mt::eTEMST_CBEI);
	cb_sim_type->setItemData(3, tr("Convergent beam electron imaging"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("ED"), mt::eTEMST_ED);
	cb_sim_type->setItemData(4, tr("Electron diffraction"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("HRTEM"), mt::eTEMST_HRTEM);
	cb_sim_type->setItemData(5, tr("High resolution transmission electron microscopy"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("PED"), mt::eTEMST_PED);
	cb_sim_type->setItemData(6, tr("Precession electron diffraction"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("HCTEM"), mt::eTEMST_HCTEM);
	cb_sim_type->setItemData(7, tr("Hollow cone transmission electron microscopy"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("EWFS"), mt::eTEMST_EWFS);
	cb_sim_type->setItemData(8, tr("Exit wave in Fourier space"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("EWRS"), mt::eTEMST_EWRS);
	cb_sim_type->setItemData(9, tr("Exit wave in real space"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("STEM_EELS"), mt::eTEMST_EELS);
	cb_sim_type->setItemData(10, tr("STEM electron energy loss spectroscopy"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("EFTEM"), mt::eTEMST_EFTEM);
	cb_sim_type->setItemData(11, tr("Energy filtered transmission electron microscopy"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("IWFS"), mt::eTEMST_IWFS);
	cb_sim_type->setItemData(12, tr("Incident wave in Fourier space"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("IWRS"), mt::eTEMST_IWRS);
	cb_sim_type->setItemData(13, tr("Incident wave in real space"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("PPFS"), mt::eTEMST_PPFS);
	cb_sim_type->setItemData(14, tr("Projected potential in the Fourier space"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("PPRS"), mt::eTEMST_PPRS);
	cb_sim_type->setItemData(15, tr("Projected potential in real space"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("TFFS"), mt::eTEMST_TFFS);
	cb_sim_type->setItemData(16, tr("Transmission function in Fourier space"), Qt::ToolTipRole);
	cb_sim_type->addItem(tr("TFRS"), mt::eTEMST_TFRS);
	cb_sim_type->setItemData(17, tr("Transmission function in real space"), Qt::ToolTipRole);
	cb_sim_type->setCurrentIndex(-1);

	disable_item_QComboBox(cb_sim_type, 12);
	disable_item_QComboBox(cb_sim_type, 14);
	disable_item_QComboBox(cb_sim_type, 16);

	auto lyh_experiment = new QHBoxLayout;
	lyh_experiment->setContentsMargins(2, 0, 1, 0);
	lyh_experiment->addWidget(lb_sim_type);
	lyh_experiment->addWidget(cb_sim_type);

	create_iw_widget();
	create_stem_widget();
	create_pcs_widget();
	create_eels_widget();
	create_pptf_widget();

	fr_experiment->layout->addLayout(lyh_experiment, 0, 0);
	fr_experiment->layout->addWidget(wg_iw, 1, 0);
	fr_experiment->layout->addWidget(wg_stem, 2, 0);
	fr_experiment->layout->addWidget(wg_pcs, 2, 0);
	fr_experiment->layout->addWidget(wg_eels, 2, 0);
	fr_experiment->layout->addWidget(wg_pptf, 2, 0);
	fr_experiment->layout->setRowStretch(3, 1);

	//////////////////////////////////////////////
	dw_simulation = new QDockWidget(tr("Experiment"), this);
	dw_simulation->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	dw_simulation->setFeatures(QDockWidget::NoDockWidgetFeatures);
	dw_simulation->setWidget(fr_experiment);
	addDockWidget(Qt::RightDockWidgetArea, dw_simulation);
}

void MainWindow::create_microscope_dock_widget()
{
	//////////////////////////////////////////////
	auto *fr_microscope = new MFrame_V;
	fr_microscope->layout->setSpacing(1);
	fr_microscope->setStyleSheet("QFrame { background: "+ QColor("#D8D8D8").name() +" }");

	/*******************************************************/
	/**********************Microscope***********************/
	/*******************************************************/
	lb_mic_name = new QLabel(tr("Microscope"));

	cb_mic_name = new QComboBox;
	cb_mic_name->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);

	cb_mic_name->setStatusTip(tr("Microscope name"));
	cb_mic_name->addItem(tr("Titan"));
	cb_mic_name->setCurrentIndex(0);

	lb_mic_acc_vol = new QLabel(tr("E<sub style=""font-family:calibri;font-size:x-large;"">0</sub> [kV]"));

	le_mic_acc_vol = new QLineEdit;
	le_mic_acc_vol->setStatusTip(tr("Acceleration voltage"));
	le_mic_acc_vol->setValidator(new QDoubleValidator(0, 5000, 1));
	le_mic_acc_vol->setText("300.0");

	pb_mic_view = new QPushButton(QIcon(":/images/view.png"), "");
	pb_mic_view->setStatusTip(tr("Add or edit microscope parameters"));
	pb_mic_view->setMaximumSize(22, 22);
	pb_mic_view->setEnabled(false);

	auto lyh_mic = new QHBoxLayout;
	lyh_mic->setSpacing(3);
	lyh_mic->setContentsMargins(2, 0, 1, 0);

	lyh_mic->addWidget(lb_mic_name);
	lyh_mic->addWidget(cb_mic_name);
	lyh_mic->addWidget(lb_mic_acc_vol);
	lyh_mic->addWidget(le_mic_acc_vol);
	lyh_mic->addWidget(pb_mic_view);

	fr_microscope->layout->addLayout(lyh_mic);

	/*******************************************************/
	/*******************Illumination model******************/
	/*******************************************************/
	lb_tilt_theta = new QLabel(tr("<span style=""font-family:symbol;font-size:large"">q</span> [°]"));
	lb_tilt_theta->setAlignment(Qt::AlignCenter);

	le_tilt_theta = new QLineEdit;
	le_tilt_theta->setStatusTip(tr("Tilt illumination polar angle"));
	le_tilt_theta->setValidator(new QDoubleValidator(0, 180, 3));
	le_tilt_theta->setText("0.000");

	/*******************************************************/
	lb_tilt_phi = new QLabel(tr("<span style=""font-family:symbol;font-size:large"">j</span> [°]"));

	le_tilt_phi = new QLineEdit;
	le_tilt_phi->setStatusTip(tr("Tilt illumination azimuthal angle"));
	le_tilt_phi->setValidator(new QDoubleValidator(0, 360, 3));
	le_tilt_phi->setText("0.000");

	/*******************************************************/
	lb_illu_model = new QLabel(tr("Illum. mode"));

	cb_illu_model = new QComboBox;
	cb_illu_model->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);

	cb_illu_model->setStatusTip(tr("Illumination mode"));
	cb_illu_model->addItem(tr("Coherence mode"), enum_to_QVar(mt::eIM_Coherent));
	cb_illu_model->setItemData(0, tr("Coherence mode illumination"), Qt::ToolTipRole);
	cb_illu_model->addItem(tr("Part. coherence mode"), enum_to_QVar(mt::eIM_Partial_Coherent));
	cb_illu_model->setItemData(1, tr("Partially coherence mode illumination"), Qt::ToolTipRole);
	cb_illu_model->addItem(tr("Trans. cross Coef."), enum_to_QVar(mt::eIM_Trans_Cross_Coef));
	cb_illu_model->setItemData(2, tr("Transmission cross coefficient"), Qt::ToolTipRole);
	cb_illu_model->addItem(tr("Incoh. numerical int."), enum_to_QVar(mt::eIM_Full_Integration));
	cb_illu_model->setItemData(3, tr("Incoherent illumination using numerical integration"), Qt::ToolTipRole);

	disable_item_QComboBox(cb_illu_model, 2, true);

	/*******************************************************/
	lb_illu_beta = new QLabel(tr("<span style=""font-family:symbol;font-size:large"">b</span> [mrad]"));
	lb_illu_beta->setAlignment(Qt::AlignCenter);

	le_illu_beta = new QLineEdit;
	le_illu_beta->setStatusTip(tr("Semi-divergence angle"));
	le_illu_beta->setValidator(new QDoubleValidator(0, std::numeric_limits<double>::max(), 3));
	le_illu_beta->setText("0.100");

	/*******************************************************/
	lb_illu_nbeta = new QLabel(tr("N<sub style=""font-family:symbol;font-size:large"">b</sub>"));
	lb_illu_nbeta->setAlignment(Qt::AlignCenter);

	le_illu_nbeta = new QLineEdit;
	le_illu_nbeta->setStatusTip(tr("Number of integration points"));
	le_illu_nbeta->setValidator(new QIntValidator(4, int_lim));
	le_illu_nbeta->setText("8");

	/*******************************************************/
	gbg_illu = new MGroupBox_G();
	gbg_illu->layout->setVerticalSpacing(1);
	gbg_illu->layout->setContentsMargins(2, 3, 1, 0);

	gbg_illu->addWidget(lb_tilt_theta, 0, 0);
	gbg_illu->addWidget(le_tilt_theta, 0, 1);
	gbg_illu->addWidget(lb_tilt_phi, 0, 2);
	gbg_illu->addWidget(le_tilt_phi, 0, 3);

	gbg_illu->addWidget(lb_illu_model, 1, 0);
	gbg_illu->addWidget(cb_illu_model, 1, 1, 1, 3);

	gbg_illu->addWidget(lb_illu_beta, 2, 0);
	gbg_illu->addWidget(le_illu_beta, 2, 1);
	gbg_illu->addWidget(lb_illu_nbeta, 2, 2);
	gbg_illu->addWidget(le_illu_nbeta, 2, 3);

	fr_microscope->layout->addWidget(gbg_illu);

	/*******************************************************/
	create_condenser_lens_tab(wg_cond_lens);
	create_objective_lens_tab(wg_obj_lens);

	tb_lens = new QTabWidget;
	tb_lens->addTab(wg_cond_lens, tr("Condenser Lens"));
	tb_lens->addTab(wg_obj_lens, tr("Objective Lens"));

	fr_microscope->layout->addWidget(tb_lens);
	fr_microscope->layout->addStretch(1);

	//////////////////////////////////////////////
	dw_microscope = new QDockWidget(tr("Microscope settings"), this);
	dw_microscope->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	dw_microscope->setFeatures(QDockWidget::NoDockWidgetFeatures);
	dw_microscope->setWidget(fr_microscope);
	addDockWidget(Qt::RightDockWidgetArea, dw_microscope);

	// signal to slot
	connect(cb_illu_model, SIGNAL(currentIndexChanged(int)), SLOT(cb_illu_model_currentIndexChanged(int)));

	cb_illu_model_currentIndexChanged(cb_illu_model->findData(mt::eIM_Coherent));
}
