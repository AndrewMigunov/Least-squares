using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;

namespace MethodKV
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        double a = -1, b = 1;
        int n = 2;
        double[] x;
        bool key = false;

        public MainWindow()
        {
            InitializeComponent();
            
        }

        private void Ok_Click(object sender, RoutedEventArgs e)
        {
            key = false;

            //this.plotView.Model = ZeroCrossing();
            //this.Ok.Content = trapezoidMethod(0, 100,function).ToString();
            //int n = 3;
            //double[,] a = new double[,] { { 1, 1, 1 }, { 2, 1, 1 }, { 1, 1, 2 } };
            //double[] y = new double[] { 1, 1, 1 };
            //double[] x = gauss(a, y, n);
            //this.Ok.Content = $"x1={x[0]},x2={x[1]},x3={x[2]}";
            double[] t = new double[n + 2];
            t[0] = a;
            t[n + 1] = b;
            double h = (b - a) / (n - 1);
            for (int j = 0; j <= n - 1; j++)
            {
                t[j + 1] = a + j * h;
            }

            double[,] matrix = new double[n, n];
            double[] y = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    matrix[i, j] = 0;
                }
            }
            matrix[0, 0] = (t[2] - t[0]) / 3;
            matrix[0, 1] = (t[2] - t[1]) / 6;
            matrix[n-1, n-2] = (t[n] - t[n-1]) / 6;
            matrix[n-1, n-1] = (t[n+1] - t[n-1]) / 3;

            for(int i = 1; i < n - 1; i++)
            {
                matrix[i, i - 1] = (t[i + 1] - t[i]) / 6;
                matrix[i, i] = (t[i + 2] - t[i]) / 3;
                matrix[i, i + 1] = (t[i + 2] - t[i + 1]) / 6;
            }

            for(int i = 0; i < n; i++)
            {
                y[i] = trapezoidMethod(a, b, (x) => (H_i(x, i + 1) * function(x)));
            }

            x = gauss(matrix, y, n);
            this.Ok.Content = $"n={n}";
            this.plotView.Model = ZeroCrossing();
            //if (key == true)
            //{
            //    if (n == 800)
            //    {
            //        return;
            //    }
            //    n += 1;                
            //    Ok_Click(sender, e);

            //}

            n += 1;


        }

        public  PlotModel ZeroCrossing()
        {
            double z = a;
            Func<double, double> fun;
            fun = function;
            var plotModel = new PlotModel();
            plotModel.PlotAreaBorderThickness = new OxyThickness(0.0);
            plotModel.PlotMargins = new OxyThickness(10);
            var LinearAxis = new LinearAxis();
            LinearAxis.Maximum = 4;
            LinearAxis.Minimum = -4;
            LinearAxis.PositionAtZeroCrossing = true;
            LinearAxis.TickStyle = TickStyle.Crossing;
            plotModel.Axes.Add(LinearAxis);
            var secondLinearAxis = new LinearAxis();
            secondLinearAxis.Maximum = 4;
            secondLinearAxis.Minimum = -4;
            secondLinearAxis.PositionAtZeroCrossing = true;
            secondLinearAxis.TickStyle = TickStyle.Crossing;
            secondLinearAxis.Position = AxisPosition.Bottom;
            plotModel.Axes.Add(secondLinearAxis);
            plotModel.Series.Add(new FunctionSeries(fun, a, b, 0.001, nameof(function)));
            var series1 = new LineSeries { Title = "Series 1" };
            double eps = 0.01;
            while (z <= b)
            {
                double y = 0;
                for(int i = 0; i < n; i++)
                {
                    y += x[i] * H_i(z, i + 1);                    
                }
                //if (Math.Abs(y - function(z)) > eps)
                //{
                //    key = true;
                //    return plotModel;

                //}
                series1.Points.Add(new DataPoint(z, y));
                z += 0.001;
            }
            
            plotModel.Series.Add(series1);

            return plotModel;
        }
        private static double function(double x)
        {
            return x*x+1;
            //return Math.Sin(10 * x);
        }
        private double H_i(double x, int i)
        {
            double[] t = new double[n + 2];
            t[0] = a;
            t[n + 1] = b;
            double h = (b - a) / (n - 1);
            for (int j = 0; j <= n - 1; j++)
            {
                t[j + 1] = a + j * h;
            }

            if (i == 1)
            {
                if ((x <= t[i + 1]) && (x >= t[i]))
                {
                    return (t[i + 1] - x) / (t[i + 1] - t[i]);
                }
                else
                {
                    return 0;
                }
            }
            if (i == n)
            {
                if ((x <= t[i]) && (x >= t[i - 1]))
                {
                    return (-t[i - 1] + x) / (t[i] - t[i - 1]);
                }
                else
                {
                    return 0;
                }
            }

            if ((x <= t[i]) && (x >= t[i - 1]))
            {
                return (-t[i - 1] + x) / (t[i] - t[i - 1]);
            }
            if ((x <= t[i + 1]) && (x >= t[i]))
            {
                return (t[i + 1] - x) / (t[i + 1] - t[i]);
            }
            else
            {
                return 0;
            }

        }

        public static double trapezoidMethod(double a, double b,Func<double,double> F)
        {
            int n = (int)((b - a) * 10000);
            double h = (b - a) / n;
            double integral = 0;
            for (int i = 0; i <= n; i++)
            {
                integral += h * F(a + i * h);
            }
            integral -= h * (F(a) + F(b));
            return integral;
        }
        public static double[] gauss(double[,] a, double[] y, int n)
        {
            double[] x;
            double max;
            int k, index;
            const double eps = 0.00001;  // точность
            x = new double[n];
            k = 0;
            while (k < n)
            {
                // Поиск строки с максимальным a[i][k]
                max = Math.Abs(a[k, k]);
                index = k;
                for (int i = k + 1; i < n; i++)
                {
                    if (Math.Abs(a[i, k]) > max)
                    {
                        max = Math.Abs(a[i, k]);
                        index = i;
                    }
                }
                // Перестановка строк
                if (max < eps)
                {
                    // нет ненулевых диагональных элементов

                    return null;
                }
                for (int j = 0; j < n; j++)
                {
                    double Temp = a[k, j];
                    a[k, j] = a[index, j];
                    a[index, j] = Temp;
                }
                double temp = y[k];
                y[k] = y[index];
                y[index] = temp;
                // Нормализация уравнений
                for (int i = k; i < n; i++)
                {
                    double Temp = a[i, k];
                    if (Math.Abs(Temp) < eps) continue; // для нулевого коэффициента пропустить
                    for (int j = 0; j < n; j++)
                        a[i, j] = a[i, j] / Temp;
                    y[i] = y[i] / Temp;
                    if (i == k) continue; // уравнение не вычитать само из себя
                    for (int j = 0; j < n; j++)
                        a[i, j] = a[i, j] - a[k, j];
                    y[i] = y[i] - y[k];
                }
                k++;
            }
            // обратная подстановка
            for (k = n - 1; k >= 0; k--)
            {
                x[k] = y[k];
                for (int i = 0; i < k; i++)
                    y[i] = y[i] - a[i, k] * x[k];
            }
            return x;
        }

        

    }
}
