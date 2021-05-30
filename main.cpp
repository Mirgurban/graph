#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>
#include <QString>

#include "window.h"

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);

  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar(window);
  Window *graph_area = new Window(window);
  QAction *action;

  if (graph_area->parse_command_line(argc, argv))
  {
    QMessageBox::warning(0, "Wrong input arguments!",
                         "Wrong input arguments!");
    return -1;
  }

  action = tool_bar->addAction("&Change function", graph_area, SLOT(change_func()));
  action->setShortcut(QString("Ctrl+C"));

  action = tool_bar->addAction("&Show 1", graph_area, SLOT(show_first_method()));
  action->setShortcut(QString("1"));

  action = tool_bar->addAction("&Show 2", graph_area, SLOT(show_second_method()));
  action->setShortcut(QString("2"));

  action = tool_bar->addAction("&Zoom +", graph_area, SLOT(zoom_in()));
  action->setShortcut(QString("+"));

  action = tool_bar->addAction("&Zoom -", graph_area, SLOT(zoom_out()));
  action->setShortcut(QString("-"));

  action = tool_bar->addAction("&+ appr points", graph_area, SLOT(enlarge_n_approx()));
  action->setShortcut(QString("Ctrl++"));

  action = tool_bar->addAction("&- appr points", graph_area, SLOT(decrise_n_approx()));
  action->setShortcut(QString("Ctrl+-"));

  action = tool_bar->addAction("&+ variation", graph_area, SLOT(enlarge_variation()));
  action->setShortcut(QString("Alt++"));

  action = tool_bar->addAction("&- variation", graph_area, SLOT(decrise_variation()));
  action->setShortcut(QString("Alt+-"));

  action = tool_bar->addAction("&Show error graph", graph_area, SLOT(show_err()));
  action->setShortcut(QString("3"));

  action = tool_bar->addAction("E&xit", window, SLOT(close()));
  action->setShortcut(QString("Ctrl+X"));

  tool_bar->setMaximumHeight(30);

  window->setMenuBar(tool_bar);
  window->setCentralWidget(graph_area);
  window->setWindowTitle("Graph");
  window->show();
  return app.exec();
}
