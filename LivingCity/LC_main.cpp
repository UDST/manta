#define BOOST_TT_HAS_OPERATOR_HPP_INCLUDED

#include <QApplication>
#include <QDebug>
#include "LC_UrbanMain.h"
#include "traffic\b18CommandLineVersion.h"

using namespace LC;

// Check command_line_options for default options.

int main(int argc, char *argv[]) {
  QSettings settings(QApplication::applicationDirPath() + "/command_line_options.ini", QSettings::IniFormat);
  bool useGUI = settings.value("GUI", false).toBool();

  if (useGUI == true) {
    QApplication a(argc, argv);
    qDebug() << "App path: " << a.applicationDirPath();
    LCUrbanMain w;
    w.showMaximized();
    return a.exec();
  } else {
    QCoreApplication a(argc, argv);
    qDebug() << "App path: " << a.applicationDirPath();

    B18CommandLineVersion cl;
    cl.runB18Simulation();
    printf(">>Simulation Ended\n");
  }
}
