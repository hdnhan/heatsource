using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace HeatSource1D_Q_fxt
{
    class Program
    {
        static void Main(string[] args)
        {
            // DATA
            Func<double, double, double> a = (x, t) => 1;

            // Ex 1:
            //Func<double, double, double> u = (x, t) => Math.Cos(Math.PI * t) * Math.Sin(Math.PI * x);
            //Func<double, double, double> F = (x, t) => -Math.PI * Math.Sin(Math.PI * t) * Math.Sin(Math.PI * x) + Math.PI * Math.PI * u(x, t);
            //Func<double, double, double> dxu = (x, t) => Math.PI * Math.Cos(Math.PI * t) * Math.Cos(Math.PI * x);

            // Ex 2:
            Func<double, double, double> u = (x, t) => Math.Exp(t) * Math.Sin(Math.PI * x);
            Func<double, double, double> F = (x, t) => (1 + Math.PI * Math.PI) * u(x, t);
            Func<double, double, double> dxu = (x, t) => Math.PI * Math.Exp(t) * Math.Cos(Math.PI * x);


            //Func<double, double, double> uD = (x, t) => 0; //have to, unuseful
            Func<double, double> u0 = (x) => u(x, 0);
            Func<double, double> dxu0 = (x) => Math.PI * Math.Cos(Math.PI * x);

            Func<double, double, double> f = (x, t) => (1 + x * x * x) * (1 + t * t);
            Func<double, double, double> q = (x, t) => 2 + t * t;
            Func<double, double, double> g = (x, t) => F(x, t) - f(x, t) * q(x, t);
            double gamma = 1e-6;


            int nn = 16;
            int Nx = nn;
            int Nt = nn;
            double[] xlim = new double[] { 0, 1 };
            double T = 1;
            double hx = (xlim[1] - xlim[0]) / Nx;
            double ht = T / Nt; // time step
            double jacobi = hx * ht; // Jacobi = 2 * area(elem)

            // Domain 2D
            double[,] node = new double[(Nx + 1) * (Nt + 1), 3];
            int[,] elem = new int[2 * Nx * Nt, 3];
            int[] dirichlet = new int[Nx + 2 * Nt + 1]; // Dirichlet boundary


            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();
            InverseFEM ifem = new InverseFEM();

            sfem.Domain(ref node, ref elem, ref dirichlet, xlim, T, Nx, Nt);


            // Observation
            double[] omega = new double[node.GetLength(0)];
            Random random = new Random();
            double eps = 0.00;
            for (int i = 0; i < node.GetLength(0); i++)
            {
                omega[i] = u(node[i, 0], node[i, 1]) - eps * (2 * random.NextDouble() - 1);
            } 

            double[] fh = ifem.CG(node, elem, dirichlet, Nx, Nt, jacobi, omega, gamma, eps, a, u0, dxu0, f, q, g);
            //double[] uh = sfem.SolveFEM(node, elem, dirichlet, jacobi, a, dxu0, oprs.Ones(node.GetLength(0)), F);

            double[] del = new double[node.GetLength(0)];
            for (int i = 0; i < node.GetLength(0); i++)
            {
                del[i] = fh[i] - f(node[i, 0], node[i, 1]);
                //del[i] =  uh[i]+ u0(node[i, 0]) - u(node[i, 0], node[i, 1]);
            }
            //double[] test = new double[] { 1, 2, 3, 4 ,5,6};
            //oprs.Print(test);
            //oprs.Print(oprs.flip(test, 1, 2));
            String path_fh = Path.Combine(Environment.CurrentDirectory, @"..\..\..\..\Results\Q_1D_fxt_" + Nx + "_" + eps + ".txt");
            File.WriteAllLines(path_fh, fh.Select(dd => dd.ToString()));
            Console.WriteLine("Error in L2: " + slvs.ErrorL2(f, fh, node, elem, jacobi));
            //Console.WriteLine("Error in W: " + slvs.ErrorW(a, dxu, uh, node, elem, jacobi));
            Console.WriteLine("min and max: " + del.Min() + ", " + del.Max());
            Console.WriteLine("Done!");
            Console.ReadLine();
        }
    }
}
