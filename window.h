
#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QMessageBox>
#include <QString>

#include "compute_polynom_in_x.h"
#include "first_method.h"
#include "second_method.h"

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  const char *f_name;
  double a;
  double b;
  int n;
  int n_appr;
  double f(double);
  double *F;
  double *X;
  double *koeff;
  double *F_2;
  double *koeff_2;
  bool show_graph_1;
  bool show_graph_2;
  bool show_graph_err;
  bool init;
  bool init_2;
  double zoom_parameter;
  double max_of_f;
  double max_of_f_1;
  double max_of_f_2;
  int variation;

public:
  Window(QWidget *parent);
  ~Window();
  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  int parse_command_line(int argc, char *argv[]);
  double Pol(double x, int n_appr);
  double Pol_2(double x, int n_appr);
  double errorfun(double x, int n_appr);
  double errorfun_2(double x, int n_appr);
  int Pol_init(int n_appr, double a, double b);
  int Pol2_init(int n_appr, double a, double b);
public slots:
  void change_func();
  void show_first_method();
  void show_second_method();
  void enlarge_n_approx();
  void decrise_n_approx();
  void enlarge_variation();
  void decrise_variation();
  void show_err();
  void zoom_in();
  void zoom_out();


protected:
  void paintEvent(QPaintEvent *event);
};

#endif
