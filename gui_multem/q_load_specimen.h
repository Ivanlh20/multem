#ifndef QT_LOAD_SPECIMEN_H
#define QT_LOAD_SPECIMEN_H

#include <QtCore>
#include "q_types.h"
#include "atom_data.hpp"

template<class T>
class Load_Specimen
{
public:
	Load_Specimen(){}
	~Load_Specimen(){}

	bool operator()(QString filename, mt::Atom_Data<double> &atoms)
	{
		QFileInfo file(filename);
		QString ext = file.suffix();

		bool result = false;
		if(ext.compare("txt", Qt::CaseInsensitive)==0)
		{
			result = read_txt(filename, atoms);
		}

		return result;
	}

private:
	bool read_txt(QString filename, mt::Atom_Data<T> &atoms)
	{
		QFile file(filename);

		if(!file.open(QFile::ReadOnly | QFile::Text)) return false;

		QTextStream text_stream(&file);

		QRegExp sep("(\\ |\\t)");
		// read first line whihc contains simulation box information
		QString line_string = text_stream.readLine().simplified();
		QStringList line_list = line_string.split(sep, QString::SkipEmptyParts);
		atoms.l_x = line_list[0].toDouble();
		atoms.l_y = line_list[1].toDouble();
		atoms.l_z = 0;
		atoms.dz = max(line_list[2].toDouble(), 0.5);

		// count number of atoms
		int n_atoms = 0;
		line_string = text_stream.readLine();
		while(!line_string.isNull())
		{
			n_atoms++;
			line_string = text_stream.readLine();
		}
		text_stream.seek(0);

		// read data
		atoms.resize(n_atoms);
		int iatoms_c = 0;
		// skip first line
		line_string = text_stream.readLine();
		// read second line
		line_string = text_stream.readLine();
		do
		{
			line_list = line_string.split(sep, QString::SkipEmptyParts);
			auto n_cols = line_list.count();

			atoms.Z[iatoms_c] = line_list[0].toDouble();
			atoms.x[iatoms_c] = line_list[1].toDouble();
			atoms.y[iatoms_c] = line_list[2].toDouble();
			atoms.z[iatoms_c] = line_list[3].toDouble();
			atoms.sigma[iatoms_c] = (n_cols>4)?line_list[4].toDouble():0.085;
			atoms.occ[iatoms_c] = (n_cols>5)?line_list[5].toDouble():1.0;
			atoms.region[iatoms_c] = (n_cols>6)?line_list[6].toDouble():0;
			atoms.charge[iatoms_c] = (n_cols>7)?line_list[7].toDouble():0;
			iatoms_c++;

			line_string = text_stream.readLine().simplified();
		} while(!line_string.isNull());

		file.close();

		atoms.get_statistic();

		return true;
	}
};

#endif // QT_LOAD_SPECIMEN_H
