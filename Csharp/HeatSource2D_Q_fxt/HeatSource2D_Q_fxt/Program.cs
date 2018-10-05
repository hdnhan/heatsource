using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_Q_fxt {
    class Program {
        static void Main(string[] args) {
            // DATA
            double axx(double x, double y, double t) => 1;
            double ayy(double x, double y, double t) => 1;

            // Ex 1:
            double u(double x, double y, double t) => Math.Cos(Math.PI * t) * Math.Sin(Math.PI * x) * Math.Sin(Math.PI * y);
            double F(double x, double y, double t) => -Math.PI * Math.Sin(Math.PI * t) * Math.Sin(Math.PI * x) * Math.Sin(Math.PI * y)
                                                                  + 2 * Math.PI * Math.PI * u(x, y, t);
            double dxu(double x, double y, double t) => Math.PI * Math.Cos(Math.PI * t) * Math.Cos(Math.PI * x) * Math.Sin(Math.PI * y);
            double dyu(double x, double y, double t) => Math.PI * Math.Cos(Math.PI * t) * Math.Sin(Math.PI * x) * Math.Cos(Math.PI * y);

            // Ex 2:


            double u0(double x, double y) => u(x, y, 0);
            double dxu0(double x, double y) => Math.PI * Math.Cos(Math.PI * x) * Math.Sin(Math.PI * y);
            double dyu0(double x, double y) => Math.PI * Math.Sin(Math.PI * x) * Math.Cos(Math.PI * y);

            double f(double x, double y, double t) => Math.Sin(Math.PI * t) * Math.Sin(Math.PI * x) * Math.Sin(Math.PI * y);
            double q(double x, double y, double t) => 1 + x * x + y * y + t * t;
            double g(double x, double y, double t) => F(x, y, t) - f(x, y, t) * q(x, y, t);
            double gamma = 1e-5;

            int nn = 2;
            int Nx = nn;
            int Ny = nn;
            int Nt = nn;
            double[] xlim = new double[] { 0, 1 };
            double[] ylim = new double[] { 0, 1 };
            double T = 1;
            double hx = (xlim[1] - xlim[0]) / Nx;
            double hy = (ylim[1] - ylim[0]) / Nt;
            double ht = T / Nt; // time step
            double jacobi = hx * hy * ht; // Jacobi = 6 * volume(elem)

            // Domain 3D
            double[,] node = new double[(Nx + 1) * (Ny + 1) * (Nt + 1), 3];
            int[,] elem = new int[6 * Nx * Ny * Nt, 4];
            int[] dirichlet = new int[(Nx + 1) * (Ny + 1) * (Nt + 1) - (Nx - 1) * (Ny - 1) * Nt]; // Dirichlet boundary

            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();
            InverseFEM ifem = new InverseFEM();
            
            //oprs.Print(omega);
            sfem.Domain(ref node, ref elem, ref dirichlet, xlim, ylim, T, Nx, Ny, Nt);

            // Observation
            double[] omega = new double[node.GetLength(0)];
            Random random = new Random();
            double eps = 0.01;
            for (int i = 0; i < node.GetLength(0); i++) {
                omega[i] = u(node[i, 0], node[i, 1], node[i, 2]) - eps * (2 * random.NextDouble() - 1);               
            }
            
            double[] fh = ifem.CG(node, elem, dirichlet, Nx, Ny, Nt, jacobi, omega, gamma, eps, axx, ayy, u0, dxu0, dyu0, f, q, g);
            double[] del = new double[node.GetLength(0)];
            for (int i = 0; i < node.GetLength(0); i++) {
                del[i] = fh[i] - f(node[i, 0], node[i, 1], node[i, 2]);              
            }
            Console.WriteLine("Error in L2: " + slvs.ErrorL2(f, fh, node, elem, jacobi));
            //Console.WriteLine("Error in W: " + slvs.ErrorW(a, dxu, uh, node, elem, jacobi));
            Console.WriteLine("min and max: " + del.Min() + ", " + del.Max());
            Console.WriteLine("Done!");
            Console.ReadLine();
        }
    }
}
