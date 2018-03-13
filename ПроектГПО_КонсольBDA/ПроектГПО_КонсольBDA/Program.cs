using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Web.UI.DataVisualization.Charting;
using System.IO;

namespace ПроектГПОКонсоль
{
    delegate double TFunc(double[] x);

    class Program
    {
        public static double[] Aver;
        public static Random rnd11 = new Random();
        static void Main(string[] args)
        {
            int SearchAgents_no, Max_iteration, dim; string Function_name;
            TFunc fobj; double Best_score = 100;
            double[] Best_p;
            SearchAgents_no = 30; // Number of search agents
            Function_name = "F1"; // Name of the test function that can be from F1 to F13 (Table 1,2,3 in the paper)
            Max_iteration = 500; // Maximum number of iterations
            int Max_progon = 30; Aver = new double[Max_iteration];
            //Load details of the selected benchmark function
            Get_Functions_details(Function_name, out dim, out fobj);
            for (int k = 0; k < Max_progon; k++)
            {
                BDA(SearchAgents_no, Max_iteration, dim, fobj, out Best_p, ref Best_score);
                Console.WriteLine("The best optimal value of the objective funciton found by BDA is : {0}", Best_score);
                //Console.WriteLine("The best solution obtained by BDA is : ");
                //for (int k1 = 0; k1 < dim; k1++)
                //{
                //    Console.WriteLine(Best_p[k1]);
                //}
            }
            using (StreamWriter sw = new StreamWriter(@"C:\Users\User\Desktop\ll.txt", false))
            {
                for (int k = 0; k < Max_iteration; k++)
                {
                    Aver[k] = Aver[k] / (double)Max_progon;
                    sw.WriteLine(Aver[k]);
                }
            }
        }
        static void Get_Functions_details(string F, out int dim, out TFunc fobj)
        {
            fobj = F1;
            dim = 75;
            switch (F)
            {
                case "F1":
                    break;
                case "F2":
                    fobj = F2;                    
                    dim = 75;
                    break;
                case "F3":
                    fobj = F3;
                    dim = 75;
                    break;
                case "F4":
                    fobj = F4;
                    dim = 75;
                    break;
                case "F5":
                    fobj = F5;                    
                    dim = 75;
                    break;
                case "F6":
                    fobj = F6;
                    dim = 75;
                    break;
                case "F7":
                    fobj = F7;
                    dim = 75;
                    break;
                case "F8":
                    fobj = F8;
                    dim = 75;
                    break;
                case "F9":
                    fobj = F9;
                    dim = 75;
                    break;
                case "F10":
                    fobj = F10;
                    dim = 75;
                    break;
                case "F11":
                    fobj = F11;
                    dim = 75;
                    break;
                case "F12":
                    fobj = F12;
                    dim = 75;
                    break;
                case "F13":
                    fobj = F13;
                    dim = 75;
                    break;
            }
        }
        static double F1(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + Math.Pow(elem, 2);
            }
            return sum;
        }
        static double F2(double[] x)
        {
            double sum = 0, prod = 1; foreach (double elem in x)
            {
                sum = sum + Math.Abs(elem);
                prod = prod * Math.Abs(elem);
            }
            return (sum + prod);
        }
        static double F3(double[] x)
        {
            int i, j, dim = x.Length; double y = 0;
            for (i = 0; i < dim; i++)
            {
                double sum = 0;
                for (j = 0; j <= i; j++)
                {
                    sum = sum + x[j];
                }
                y = y + Math.Pow(sum, 2);
            }
            return y;
        }
        static double F4(double[] x)
        {
            double max = 0;
            foreach (double elem in x)
            {
                if (max < Math.Abs(elem)) max = Math.Abs(elem);
            }
            return max;
        }
        static double F5(double[] x)
        {
            int dim = x.Length, i; double sum = 0;
            for (i = 0; i <= (dim - 2); i++)
            {
                sum = sum + (100 * Math.Pow((x[i + 1] - Math.Pow(x[i], 2)), 2) + Math.Pow((x[i] - 1), 2));
            }
            return sum;
        }
        static double F6(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + Math.Pow((elem + 0.5), 2);
            }
            return sum;
        }
        static double F7(double[] x)
        {
            double sum = rnd11.NextDouble(); int i, dim = x.Length;
            for (i = 0; i < dim; i++)
            {
                sum = sum + (i + 1) * Math.Pow(x[i], 4);
            }
            return sum;
        }
        static double F8(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + (-1 * elem * Math.Sin(Math.Sqrt(Math.Abs(elem))));
            }
            sum = sum / (double)x.Length;
            return sum;
        }
        static double F9(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + (Math.Pow(elem, 2) - 10 * Math.Cos(2 * Math.PI * elem) + 10);
            }
            return sum;
        }
        static double F10(double[] x)
        {
            double sum1 = 0, sum2 = 0, sum;
            foreach (double elem in x)
            {
                sum1 = sum1 + Math.Pow(elem, 2);
                sum2 = sum2 + Math.Cos(2 * Math.PI * elem);
            }
            sum1 = sum1 / x.Length; sum2 = sum2 / x.Length;
            sum = -20 * Math.Exp(-0.2 * Math.Sqrt(sum1)) - Math.Exp(sum2) + 20 + Math.Exp(1);
            return sum;
        }
        static double F11(double[] x)
        {
            double sum = 0, mply = 1; int i, dim = x.Length;
            for (i = 0; i < dim; i++)
            {
                sum = sum + Math.Pow(x[i], 2);
                mply = mply * Math.Cos(x[i] / Math.Sqrt(i + 1));
            }
            double y = sum / 4000 - mply + 1;
            return y;
        }
        static double F12(double[] x)
        {
            int dim = x.Length;
            double[] z = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                z[i] = 1 + (x[i] + 1) / 4.0;
            }
            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < (dim - 1); i++)
            {
                sum1 = sum1 + Math.Pow(z[i] - 1, 2) * (1 + 10 * Math.Pow(Math.Sin(Math.PI * z[i + 1]), 2));
            }
            foreach (double elem in x)
            {
                sum2 = sum2 + uFunc(elem, 10, 100, 4);
            }
            double y = (Math.PI / dim) * (10 * Math.Sin(Math.PI * z[0]) + sum1 + Math.Pow(z[dim - 1] - 1, 2)) + sum2;
            return y;
        }
        static double uFunc(double x, int a, int k, int m)
        {
            if (x > a) return (k * Math.Pow(x - a, m));
            if (x < -a) return (k * Math.Pow(-x - a, m));
            return 0;
        }
        static double F13(double[] x)
        {
            int dim = x.Length; double sum1 = 0, sum2 = 0;
            foreach (double elem in x)
            {
                sum1 = sum1 + Math.Pow(elem - 1, 2) * (1 + Math.Pow(Math.Sin(3 * Math.PI * elem + 1), 2));
                sum2 = sum2 + uFunc(elem, 5, 100, 4);
            }
            double y = 0.1 * (Math.Pow(Math.Sin(3 * Math.PI * x[0]), 2) + sum1 + Math.Pow(x[dim - 1] - 1, 2) * (1 + Math.Pow(Math.Sin(2 * Math.PI * x[dim - 1]), 2))) + sum2;
            return y;
        }
        static double Erf(double x)
        {
            // constants
            double a1 = 0.254829592;
            double a2 = -0.284496736;
            double a3 = 1.421413741;
            double a4 = -1.453152027;
            double a5 = 1.061405429;
            double p = 0.3275911;
            // Save the sign of x
            int sign = 1;
            if (x < 0)
                sign = -1;
            x = Math.Abs(x);
            // A&S formula 7.1.26
            double t = 1.0 / (1.0 + p * x);
            double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);
            return sign * y;
        }
        static double Sigmoid(double x, double k)
        {
            return 1 / (1 + Math.Exp(-k*x));
        }
        static void initialization(int SearchAgents_no, int dim, ref double[][] dX)
        {
            int i, j;
            for (i = 0; i < SearchAgents_no; i++)
            {
                for (j = 0; j < dim; j++)
                {
                    double r1 = rnd11.NextDouble();                    
                    if (r1 <= 0.5)
                    dX[i][j] = 0;
                    else dX[i][j] = 1;
                }
            }
        }
        static void BDA(int SearchAgents_no, int Max_iteration, int dim, TFunc fobj, out double[] Best_p, ref double Best_s)
        {
            int i, j; Best_p = new double[dim];
            double Food_fitness, Enemy_fitness;
            double[] Food_pos = new double[dim], Enemy_pos = new double[dim];
            Food_fitness = Double.PositiveInfinity; Enemy_fitness = Double.NegativeInfinity;
            double[][] X = new double[SearchAgents_no][]; double[][] DeltaX = new double[SearchAgents_no][];
            for (i = 0; i < SearchAgents_no; i++)
            {
                X[i] = new double[dim];
                DeltaX[i] = new double[dim];
            }
            initialization(SearchAgents_no, dim, ref X);            
            double[] Fitness = new double[SearchAgents_no];
            initialization(SearchAgents_no, dim, ref DeltaX);            
            int iter; double my_c, w;
            for (iter = 1; iter <= Max_iteration; iter++)
            {
                w = 0.9 - iter * (0.5 / (double)Max_iteration);
                my_c = 0.1 - iter * (0.1 / ((double)Max_iteration / 2.0));
                if (my_c < 0) my_c = 0;
                double s, a, c, f, e;
                s = 2 * rnd11.NextDouble() * my_c; // Seperation weight                
                a = 2 * rnd11.NextDouble() * my_c; // Alignment weight                
                c = 2 * rnd11.NextDouble() * my_c; // Cohesion weight         
                f = 2 * rnd11.NextDouble();        // Food attraction weight
                e = my_c;                          // Enemy distraction weight
                if (iter > (3.0 * Max_iteration / 4.0))
                {
                    e = 0;
                }
                for (i = 0; i < SearchAgents_no; i++) //Calculate all the objective values
                {
                    Fitness[i] = fobj(X[i]);
                    if (Fitness[i] < Food_fitness)
                    {
                        Food_fitness = Fitness[i];
                        for (j = 0; j < dim; j++)
                        {
                            Food_pos[j] = X[i][j];
                        }
                    }
                    if (Fitness[i] > Enemy_fitness)
                    {
                            Enemy_fitness = Fitness[i];
                            for (j = 0; j < dim; j++)
                            {
                                Enemy_pos[j] = X[i][j];
                            }
                    }
                }
                for (i = 0; i < SearchAgents_no; i++)
                {
                    int index = 0; int neighbours_no = 0;
                    double[][] Neighbours_DeltaX = new double[SearchAgents_no][]; //clear?
                    double[][] Neighbours_X = new double[SearchAgents_no][];
                    //Find the neighbouring solutions (all the dragonflies are assumed as a group in binary seach spaces)
                    for (j = 0; j < SearchAgents_no; j++)
                    {
                        if (i != j)
                        {
                            index = index + 1;
                            neighbours_no = neighbours_no + 1;
                            Neighbours_DeltaX[index - 1] = new double[dim];
                            Neighbours_X[index - 1] = new double[dim];
                            for (int k = 0; k < dim; k++)
                            {
                                Neighbours_DeltaX[index - 1][k] = DeltaX[j][k];
                                Neighbours_X[index - 1][k] = X[j][k];
                            }
                        }
                    }
                    // Seperation Eq. (2.1)
                    double[] S = new double[dim];                    
                        for (int k = 0; k < neighbours_no; k++)
                        {
                            for (int ka = 0; ka < dim; ka++)
                            {
                                S[ka] = S[ka] + (Neighbours_X[k][ka] - X[i][ka]);  // mistakable?
                            }
                        }
                        for (int k = 0; k < dim; k++)
                        {
                            S[k] = -1 * S[k];
                        }
                    
                    // Alignment Eq. (2.2)
                    double[] A = new double[dim];                    
                        for (int k1 = 0; k1 < dim; k1++)
                        {
                            for (int k2 = 0; k2 < neighbours_no; k2++)
                            {
                                A[k1] = A[k1] + Neighbours_DeltaX[k2][k1];
                            }
                            A[k1] = A[k1] / (double)neighbours_no;
                        }
                    
                    // Cohesion Eq. (2.3)
                    double[] C_temp = new double[dim];                    
                        for (int k1 = 0; k1 < dim; k1++)
                        {
                            for (int k2 = 0; k2 < neighbours_no; k2++)
                            {
                                C_temp[k1] = C_temp[k1] + Neighbours_X[k2][k1];
                            }
                            C_temp[k1] = C_temp[k1] / (double)neighbours_no;
                        }                    
                    double[] C = new double[dim];
                    for (int k = 0; k < dim; k++)
                    {
                        C[k] = C_temp[k] - X[i][k];
                    }

                    // Attraction to food Eq. (2.4)                    
                    double[] F = new double[dim];            
                        for (int k = 0; k < dim; k++)
                        {
                            F[k] = Food_pos[k] - X[i][k];
                        }
                    
                    // Distraction from enemy Eq. (3.5)                    
                    double[] Enemy = new double[dim];                    
                        for (int k = 0; k < dim; k++)
                        {
                            Enemy[k] = Enemy_pos[k] + X[i][k];
                        }
                    
                    for (j = 0; j < dim; j++)
                        {
                            // Eq. (2.6)
                            DeltaX[i][j] = (s * S[j] + a * A[j] + c * C[j] + f * F[j] + e * Enemy[j]) + w * DeltaX[i][j];
                            if (DeltaX[i][j] > 6)
                            {
                                DeltaX[i][j] = 6;
                            }
                            if (DeltaX[i][j] < -6)
                            {
                                DeltaX[i][j] = -6;
                            }
                            //Transfer function  Eq. (2.11)
                            double T,r2;
                        T = Math.Abs(Erf((Math.Sqrt(Math.PI) / 2.0) * DeltaX[i][j]));
                        //T = Math.Abs(Math.Tanh(DeltaX[i][j]));
                        //T = Math.Abs(DeltaX[i][j] / Math.Sqrt(1 + Math.Pow(DeltaX[i][j], 2)));
                        //T = Math.Abs((2.0/Math.PI)*Math.Atan((Math.PI/2.0)*DeltaX[i][j]));
                        //T = Sigmoid(DeltaX[i][j], 1.0);
                        //Eq. (3.12)
                            r2 = rnd11.NextDouble();
                        if (r2 < T)
                            {
                                if (X[i][j] == 0) X[i][j] = 1;
                                else X[i][j] = 0;
                            }
                        //if (r2 < T) X[i][j] = 0; else X[i][j] = 1;
                        }
                    }
                Best_s = Food_fitness;
                Aver[iter - 1] = Aver[iter - 1] + Food_fitness;
                Best_p = Food_pos;
            }
        }
    }
}