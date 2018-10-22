using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_1PointOnOmega_ft {
    class Program {
        static void Main(string[] args) {
            // DATA
            double axx(double x, double y, double t) => 1;
            double ayy(double x, double y, double t) => 1;

            // Example:
            double u(double x, double y, double t) => Math.Exp(t) * x * (1 - x) * Math.Sin(Math.PI * y);
            double F(double x, double y, double t) => u(x, y, t) + 2 * Math.Exp(t) * Math.Sin(Math.PI * y) + Math.PI * Math.PI * u(x, y, t);

            double u0(double x, double y) => u(x, y, 0);
            double dxu0(double x, double y) => (1 - 2 * x) * Math.Sin(Math.PI * y);
            double dyu0(double x, double y) => Math.PI * x * (1 - x) * Math.Cos(Math.PI * y);

            //double f(double t) => Math.Sin(Math.PI * t);
            //double f(double t) => 2 * t * ((0.5 - t) >= 0 ? 1 : 0) + 2 * (1 - t) * ((t - 0.5) > 0 ? 1 : 0);
            double f(double t) => ((t - 0.25) >= 0 ? 1 : 0) * ((t - 0.75) <= 0 ? 1 : 0);

            double q(double x, double y, double t) => x * y + t + 1;
            double g(double x, double y, double t) => F(x, y, t) - f(t) * q(x, y, t);
            double gamma = 1e-5;

            int nn = 64;
            int Nx = nn;
            int Ny = nn;
            int Nt = nn;
            double[] xlim = new double[] { 0, 1 };
            double[] ylim = new double[] { 0, 1 };
            double T = 1;
            double hx = (xlim[1] - xlim[0]) / Nx;
            double hy = (ylim[1] - ylim[0]) / Ny;
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

            // Integral observation
            double x0 = 0.48;
            double y0 = 0.48;
            double w(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0)) * Math.Cosh(Nx * (x - x0))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0)) * Math.Cosh(Ny * (y - y0))); // w(x, y) only
            double[] omega = new double[Nt + 1];
            Random random = new Random();
            // avoid initial condition noise
            for (int nt = 1; nt <= Nt; nt++) {
                omega[nt] = 2 * random.NextDouble() - 1;
            }
            double norm_noise = Math.Sqrt(slvs.NormL2(omega, ht));
            double eps = 0.01;
            for (int nt = 0; nt <= Nt; nt++) {
                omega[nt] = u(x0, y0, ht * nt) - eps * omega[nt] / norm_noise;
            }

            double[] fh = ifem.CG(node, elem, dirichlet, Nx, Ny, Nt, T, hx, hy, ht, jacobi, omega, gamma, eps, w, axx, ayy, u0, dxu0, dyu0, f, q, g);
            //Console.ReadLine();
        }
    }
}
