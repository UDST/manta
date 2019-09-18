#define BOOST_TT_HAS_OPERATOR_HPP_INCLUDED


#ifdef B18_RUN_WITH_GUI
#include <QApplication>
#include "LC_UrbanMain.h"
#else
#include "qcoreapplication.h"
#endif
#include <QDebug>
#include "traffic/b18CommandLineVersion.h"

int main(int argc, char *argv[]) {
  QCoreApplication a(argc, argv);
  QSettings settings(QCoreApplication::applicationDirPath() +
                     "/command_line_options.ini", QSettings::IniFormat);

  LC::B18CommandLineVersion cl;
  cl.runB18Simulation();

  return 0;
}
