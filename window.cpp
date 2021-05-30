
#include <QPainter>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "window.h"

#define DEFAULT_A -10
#define DEFAULT_B 10

static double f_0(double x)
{
  return 1;
}

static double f_1(double x)
{
  return x;
}

static double f_2(double x)
{
  return x * x;
}

static double f_3(double x)
{
  return x * x * x;
}

static double f_4(double x)
{
  return x * x * x * x;
}

static double f_5(double x)
{
  return exp(x);
}

static double f_6(double x)
{
  return 1 / (25 * x * x + 1);
}


double Window::f(double x)
{

  switch (func_id)
  {
  case 0:
    f_name = "k = 0   f (x) = 1";
    if (fabs(x - a - ((b - a) / (n_appr - 1)) * (n_appr / 2)) < 1e-13)
      return f_0(x) + variation * 0.1 * max_of_f;
    return (f_0(x));
  case 1:
    f_name = "k = 1   f (x) = x";
    if (fabs(x - a - ((b - a) / (n_appr - 1)) * (n_appr / 2)) < 1e-13)
      return f_1(x) + variation * 0.1 * max_of_f;
    return (f_1(x));
  case 2:
    f_name = "k = 2   f (x) = x * x";
    if (fabs(x - a - ((b - a) / (n_appr - 1)) * (n_appr / 2)) < 1e-13)
      return f_2(x) + variation * 0.1 * max_of_f;
    return (f_2(x));
  case 3:
    f_name = "k = 3   f (x) = x * x * x";
    if (fabs(x - a - ((b - a) / (n_appr - 1)) * (n_appr / 2)) < 1e-13)
      return f_3(x) + variation * 0.1 * max_of_f;
    return (f_3(x));
  case 4:
    f_name = "k = 4   f (x) = x * x * x * x";
    if (fabs(x - a - ((b - a) / (n_appr - 1)) * (n_appr / 2)) < 1e-13)
      return f_4(x) + variation * 0.1 * max_of_f;
    return (f_4(x));
  case 5:
    f_name = "k = 5   f (x) = exp(x)";
    if (fabs(x - a - ((b - a) / (n_appr - 1)) * (n_appr / 2)) < 1e-13)
      return f_5(x) + variation * 0.1 * max_of_f;
    return (f_5(x));
  case 6:
    f_name = "k = 6   f (x) = 1 / (25 * x * x + 1)";
    if (fabs(x - a - ((b - a) / (n_appr - 1)) * (n_appr / 2)) < 1e-13)
      return f_6(x) + variation * 0.1 * max_of_f;
    return (f_6(x));
  }
  return 0;
}

int Window::Pol_init(int n_appr, double a, double b)
{
  double h = (b - a) / (double)(n_appr - 1);
  F = (double *)realloc(F, n_appr * sizeof(double));
  X = (double *)realloc(X, n_appr * sizeof(double));
  double derivative1;
  double derivativeN;
  koeff = (double *)realloc(koeff, 4 * n_appr * sizeof(double));
  for (int i = 0; i < 4 * n_appr; i++)
  {
    koeff[i] = 0;
  }
  for (int i = 0; i < n_appr; i++)
  {
    X[i] = a + i * h;
    F[i] = f(X[i]);
  }
  derivative1 = (f(a + h / 1e10) - f(a)) / (h / 1e10);
  derivativeN = (f(b) - f(b - h / 1e10)) / (h / 1e10);
  first_method(n_appr, a, b, X, F, koeff, derivative1, derivativeN);
  return 0;
}

int Window::Pol2_init(int n_appr, double a, double b)
{
    double h = (b - a) / (double)(n_appr - 1);
    F = (double *)realloc(F, n_appr * sizeof(double));
    X = (double *)realloc(X, n_appr * sizeof(double));
    double derivative1;
    double derivativeN;
    koeff_2 = (double *)realloc(koeff_2, 4 * n_appr * sizeof(double));
    for (int i = 0; i < 4 * n_appr; i++)
    {
      koeff_2[i] = 0;
    }
    for (int i = 0; i < n_appr; i++)
    {
      X[i] = a + i * h;
      F[i] = f(X[i]);
    }
    derivative1 = (f(a + h / 1e10) - f(a)) / (h / 1e10);
    derivativeN = (f(b) - f(b - h / 1e10)) / (h / 1e10);
    second_method(n_appr, X, F, koeff_2, derivative1, derivativeN);
    return 0;
}

double Window::errorfun(double x, int n_appr)
{
  return fabs(Pol(x, n_appr) - f(x));
}

double Window::errorfun_2(double x, int n_appr)
{
  return fabs(Pol_2(x, n_appr) - f(x));
}

double Window::Pol(double x, int n_appr)
{
  return compute_pol(n_appr, a, b, x, koeff);
}

double Window::Pol_2(double x, int n_appr)
{
  return compute_pol(n_appr, a, b, x, koeff_2);
}


Window::Window(QWidget *parent)
    : QWidget(parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = width();
  n_appr = 20;
  show_graph_1 = 1;
  show_graph_2 = 0;
  show_graph_err = 0;
  zoom_parameter = 1;
  variation = 0;
  init = 0;
  init_2 = 0;
  max_of_f = 0;
  max_of_f_1 = 0;
  max_of_f_2 = 0;
  func_id = 6;

  change_func();

  X = (double *)malloc(n_appr * sizeof(double));
  F = (double *)malloc(n_appr * sizeof(double));
  koeff = (double *)malloc(4 * n_appr * sizeof(double));
  F_2 = (double *)malloc(n_appr * sizeof(double));
  koeff_2 = (double *)malloc(4 * n_appr * sizeof(double));
}

Window::~Window()
{
  free(X);
  free(F);
  free(koeff);
  free(F_2);
  free(koeff_2);
}
QSize Window::minimumSizeHint() const
{
  return QSize(100, 100);
}

QSize Window::sizeHint() const
{
  return QSize(1000, 1000);
}

int Window::parse_command_line(int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
    return -1;

  if (sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 || sscanf(argv[4], "%d", &func_id) != 1 || b - a < 1.e-6 || (sscanf(argv[3], "%d", &n_appr) != 1) || n_appr <= 9 || func_id < 0 || func_id > 6)
    return -2;
  return 0;
}

void Window::change_func()
{
  func_id = (func_id + 1) % 7;

  switch (func_id)
  {
  case 0:
    f_name = "f (x) = 1";
    break;
  case 1:
    f_name = "f (x) = x";
    break;
  case 2:
    f_name = "f (x) = x * x";
    break;
  case 3:
    f_name = "f (x) = x * x * x";
    break;
  case 4:
    f_name = "f (x) = x * x * x * x";
    break;
  case 5:
    f_name = "f (x) = exp(x)";
    break;
  case 6:
    f_name = "f (x) = 1 / (25 * x * x + 1)";
    break;
  }
  init = 0;
  init_2 = 0;
  update();
}


void Window::enlarge_variation()
{
  variation += 1;
  init = 0;
  init_2 = 0;
  update();
}

void Window::decrise_variation()
{
  variation -= 1;
  init = 0;
  init_2 = 0;
  update();
}


void Window::enlarge_n_approx()
{
  n_appr += 10;
  init = 0;
  init_2 = 0;
  update();
}

void Window::decrise_n_approx()
{
  if (n_appr - 10 > 9)
  {
    n_appr -= 10;
  }
  init = 0;
  init_2 = 0;
  update();
}


void Window::show_first_method()
{
  show_graph_1 = (show_graph_1 + 1) % 2;
  init = 0;
  update();
}
void Window::show_second_method()
{
  show_graph_2 = (show_graph_2 + 1) % 2;
  init_2 = 0;
  update();
}
void Window::show_err()
{
  show_graph_err = (show_graph_err + 1) % 2;
  show_graph_1 = 0;
  show_graph_2 = 0;
  init = 0;
  init_2 = 0;
  update();
}

void Window::zoom_in()
{
  zoom_parameter *= 2;
  update();
}
void Window::zoom_out()
{
  if (zoom_parameter > 1)
  {
    zoom_parameter /= 2;
    update();
  }
}

void Window::paintEvent(QPaintEvent * )
{
  QPainter painter(this);
  double x1, x2, y1, y2;
  double max_y, min_y;
  n = width();
  double delta_y, delta_x = (b - a) / (n * zoom_parameter);


  max_y = min_y = 0;
  for (x1 = a; x1 - b < 1.e-6; x1 += delta_x)
  {
    y1 = f(x1);
    if (y1 < min_y)
      min_y = y1;
    if (y1 > max_y)
      max_y = y1;
  }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;
  painter.setPen(QPen(Qt::black, 10/height()));
  painter.save();

  painter.translate(0.5 * width(), 0.5 * height());
  painter.scale(zoom_parameter * width() / (b - a), -height() / (max_y - min_y));
  painter.translate(-0.5 * (a + b), -0.5 * (min_y + max_y));


  if (show_graph_err == 0)
  {
    max_of_f = 0;
    x1 = a;
    y1 = f(x1);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
    {
      y2 = f(x2);
      painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
      if (max_of_f < fabs(y2))
      {
        max_of_f = fabs(y2);
      }
      x1 = x2, y1 = y2;
    }
    x2 = b;
    y2 = f(x2);
    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
  }

  if (show_graph_err == 1)
  {
    max_of_f_1 = 0;
    painter.setPen(QPen(Qt::blue, 10/height()));
    if (init == 0)
    {
      Pol_init(n_appr, a, b);
      init = 1;
    }

    x1 = a;
    y1 = errorfun(x1, n_appr);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
    {
      y2 = errorfun(x2, n_appr);
      if (max_of_f_1 < y2)
      {
        max_of_f_1 = y2;
      }
      painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

      x1 = x2, y1 = y2;
    }
    x2 = b;
    y2 = errorfun(x2, n_appr);
    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

    max_of_f_2 = 0;
    painter.setPen(QPen(Qt::darkGreen, 10/height()));
    if (init_2 == 0)
    {
      Pol2_init(n_appr, a, b);
      init_2 = 1;
    }
    x1 = a;
    y1 = errorfun_2(x1, n_appr);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
    {
      y2 = errorfun_2(x2, n_appr);
      if (max_of_f_2 < y2)
      {
        max_of_f_2 = y2;
      }
      painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

      x1 = x2, y1 = y2;
    }
    x2 = b;
    y2 = errorfun_2(x2, n_appr);
    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
  }

  if (show_graph_1 == 1)
  {
    painter.setPen(QPen(Qt::blue, 10/height()));
    if (init == 0)
    {
      Pol_init(n_appr, a, b);
      init = 1;
    }
    x1 = a;
    y1 = Pol(x1, n_appr);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
    {
      y2 = Pol(x2, n_appr);
      painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

      x1 = x2, y1 = y2;
    }
//    x2 = b;
//    y2 = Pol(x2, n_appr);
//    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
  }

  if (show_graph_2 == 1)
  {
    painter.setPen(QPen(Qt::darkMagenta, 10/height()));
    if (init_2 == 0)
    {
      Pol2_init(n_appr, a, b);
      init_2 = 1;
    }
    x1 = a;
    y1 = Pol_2(x1, n_appr);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x)
    {
      y2 = Pol_2(x2, n_appr);
      painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

      x1 = x2, y1 = y2;
    }
//    x2 = b;
//    y2 = Pol_2(x2, n_appr);
//    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
  }



  painter.setPen(QPen(Qt::red, 10/height()));
  painter.drawLine(a, 0, b, 0);
  painter.drawLine(0, max_y, 0, min_y);

  painter.restore();

  painter.setPen(QPen(Qt::black, 10/height()));
  painter.drawText(0, 20, f_name);
  painter.drawText(0, 40, "approx points = " + QString::number(n_appr));

  if (show_graph_1 == 1 || show_graph_err == 1)
  {
    painter.setPen(QPen(Qt::blue, 10/height()));
    painter.drawText(0, 80, "- method 1");
    if (show_graph_err == 1)
    {
      painter.drawText(0, 100, "Err1 = " + QString::number(max_of_f_1));
    }
  }
  if (show_graph_2 == 1 || show_graph_err == 1)
  {
    painter.setPen(QPen(Qt::darkMagenta, 10/height()));
    painter.drawText(0, 140, "- method 2");
    if (show_graph_err == 1)
    {
      painter.drawText(0, 160, "Err2 = " + QString::number(max_of_f_2));
    }
  }

  painter.drawText(0, 200, "variation = " + QString::number(variation));
}
