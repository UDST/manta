#define BOOST_TT_HAS_OPERATOR_HPP_INCLUDED

#ifdef B18_RUN_WITH_GUI
#include <QApplication>
#include "LC_UrbanMain.h"
#else
#include "qcoreapplication.h"
#endif
#include <QDebug>
#include "traffic/b18CommandLineVersion.h"

// NOTE: Check command_line_options for default options.

int main(int argc, char *argv[]) {
#ifdef B18_RUN_WITH_GUI
  QApplication a(argc, argv);
  QSettings settings(QApplication::applicationDirPath() +
                     "/command_line_options.ini", QSettings::IniFormat);
  bool useGUI = settings.value("GUI", true).toBool();

  if (useGUI == true) {
    qDebug() << "App path: " << a.applicationDirPath();
    LC::LCUrbanMain w;
    w.showMaximized();
    return a.exec();
  } else {
    qDebug() << "App path: " << a.applicationDirPath();

    LC::B18CommandLineVersion cl;
    cl.runB18Simulation();
    printf(">>Simulation Ended\n");
  }

#else
  QCoreApplication a(argc, argv);
  QSettings settings(QCoreApplication::applicationDirPath() +
                     "/command_line_options.ini", QSettings::IniFormat);

  qDebug() << "App path: " << a.applicationDirPath();

  LC::B18CommandLineVersion cl;
  cl.runB18Simulation();
  printf(">>Simulation Ended\n");
#endif

}
