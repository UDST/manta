#define BOOST_TT_HAS_OPERATOR_HPP_INCLUDED

#include <QApplication>
#include <QDebug>
#include "LC_UrbanMain.h"
#include "traffic\b18CommandLineVersion.h"

using namespace LC;

const bool kUseGUI = true; // flag to control -> Run command line (false) or GUI version (true).

int main(int argc, char *argv[]) {
  QApplication a(argc, argv);
  qDebug() << "App path: " << a.applicationDirPath();
  if (kUseGUI == true) {
    LCUrbanMain w;
    w.showMaximized();
    return a.exec();
  } else {
    B18CommandLineVersion cl;
    cl.runB18Simulation();
    printf(">>Simulation Ended\n");
  }
}
