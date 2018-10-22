using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_9PointOnOmega_ft {
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

            double f(double x, double y) => Math.Sin(Math.PI * x) * y * (1 - y);
            double q(double x, double y, double t) => x * y + t + 1;
            double g(double x, double y, double t) => F(x, y, t) - f(x, y) * q(x, y, t);
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
            double[] x0 = new double[] { 0.1, 0.2, 0.4, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
            double[] y0 = new double[] { 0.4, 0.7, 0.1, 0.4, 0.9, 0.6, 0.2, 0.8, 0.5 };
            double w0(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[0])) * Math.Cosh(Nx * (x - x0[0]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[0])) * Math.Cosh(Ny * (y - y0[0]))); // w(x, y) only
            double w1(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[1])) * Math.Cosh(Nx * (x - x0[1]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[1])) * Math.Cosh(Ny * (y - y0[1])));
            double w2(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[2])) * Math.Cosh(Nx * (x - x0[2]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[2])) * Math.Cosh(Ny * (y - y0[2])));
            double w3(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[3])) * Math.Cosh(Nx * (x - x0[3]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[3])) * Math.Cosh(Ny * (y - y0[3])));
            double w4(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[4])) * Math.Cosh(Nx * (x - x0[4]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[4])) * Math.Cosh(Ny * (y - y0[4])));
            double w5(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[5])) * Math.Cosh(Nx * (x - x0[5]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[5])) * Math.Cosh(Ny * (y - y0[5])));
            double w6(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[6])) * Math.Cosh(Nx * (x - x0[6]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[6])) * Math.Cosh(Ny * (y - y0[6])));
            double w7(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[7])) * Math.Cosh(Nx * (x - x0[7]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[7])) * Math.Cosh(Ny * (y - y0[7])));
            double w8(double x, double y, double t) => Nx / (2 * Math.Cosh(Nx * (x - x0[8])) * Math.Cosh(Nx * (x - x0[8]))) *
                                                      Ny / (2 * Math.Cosh(Ny * (y - y0[8])) * Math.Cosh(Ny * (y - y0[8])));



            double[] omega0 = new double[Nt + 1];
            double[] omega1 = new double[Nt + 1];
            double[] omega2 = new double[Nt + 1];
            double[] omega3 = new double[Nt + 1];
            double[] omega4 = new double[Nt + 1];
            double[] omega5 = new double[Nt + 1];
            double[] omega6 = new double[Nt + 1];
            double[] omega7 = new double[Nt + 1];
            double[] omega8 = new double[Nt + 1];
            Random random = new Random();
            // avoid initial condition noise

            for (int nt = 1; nt <= Nt; nt++) {
                omega0[nt] = 2 * random.NextDouble() - 1;
                omega1[nt] = 2 * random.NextDouble() - 1;
                omega2[nt] = 2 * random.NextDouble() - 1;
                omega3[nt] = 2 * random.NextDouble() - 1;
                omega4[nt] = 2 * random.NextDouble() - 1;
                omega5[nt] = 2 * random.NextDouble() - 1;
                omega6[nt] = 2 * random.NextDouble() - 1;
                omega7[nt] = 2 * random.NextDouble() - 1;
                omega8[nt] = 2 * random.NextDouble() - 1;
            }
            double[] norm_noise = new double[] { Math.Sqrt(slvs.NormL2(omega0, ht)),
                                                  Math.Sqrt(slvs.NormL2(omega1, ht)),
                                                  Math.Sqrt(slvs.NormL2(omega2, ht)),
                                                  Math.Sqrt(slvs.NormL2(omega3, ht)),
                                                  Math.Sqrt(slvs.NormL2(omega4, ht)),
                                                  Math.Sqrt(slvs.NormL2(omega5, ht)),
                                                  Math.Sqrt(slvs.NormL2(omega6, ht)),
                                                  Math.Sqrt(slvs.NormL2(omega7, ht)),
                                                  Math.Sqrt(slvs.NormL2(omega8, ht))};

            double eps = 0.01;
            
            for (int nt = 0; nt <= Nt; nt++) {
                omega0[nt] = u(x0[0], y0[0], ht * nt) - eps * omega0[nt] / norm_noise[0];
                omega1[nt] = u(x0[1], y0[1], ht * nt) - eps * omega1[nt] / norm_noise[1];
                omega2[nt] = u(x0[2], y0[2], ht * nt) - eps * omega2[nt] / norm_noise[2];
                omega3[nt] = u(x0[3], y0[3], ht * nt) - eps * omega3[nt] / norm_noise[3];
                omega4[nt] = u(x0[4], y0[4], ht * nt) - eps * omega4[nt] / norm_noise[4];
                omega5[nt] = u(x0[5], y0[5], ht * nt) - eps * omega5[nt] / norm_noise[5];
                omega6[nt] = u(x0[6], y0[6], ht * nt) - eps * omega6[nt] / norm_noise[6];
                omega7[nt] = u(x0[7], y0[7], ht * nt) - eps * omega7[nt] / norm_noise[7];
                omega8[nt] = u(x0[8], y0[8], ht * nt) - eps * omega8[nt] / norm_noise[8];
            }

            double[] fh = ifem.CG(node, elem, dirichlet, Nx, Ny, Nt, T, hx, hy, ht, jacobi, gamma, eps, 
                omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7, omega8, 
                w0, w1, w2, w3, w4, w5, w6, w7, w8, axx, ayy, u0, dxu0, dyu0, f, q, g);
            //Console.ReadLine();
        }
    }
}
