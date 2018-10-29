using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_IntegrationOnOmega_ft {
    class Program {
        static void Main(string[] args) {
            // DATA
            double axx(double x, double y, double t) => 1;
            double ayy(double x, double y, double t) => 1;

            // Ex 1:
            double u(double x, double y, double t) => Math.Exp(t) * x * (1 - x) * Math.Sin(Math.PI * y);
            double F(double x, double y, double t) => u(x, y, t) + 2 * Math.Exp(t) * Math.Sin(Math.PI * y) + Math.PI * Math.PI * u(x, y, t);

            double u0(double x, double y) => u(x, y, 0);
            double dxu0(double x, double y) => (1 - 2 * x) * Math.Sin(Math.PI * y);
            double dyu0(double x, double y) => Math.PI * x * (1 - x) * Math.Cos(Math.PI * y);

            //double f(double t) => Math.Sin(Math.PI * t);
            double f(double t) => 2 * t * ((0.5 - t) >= 0 ? 1 : 0) + 2 * (1 - t) * ((t - 0.5) > 0 ? 1 : 0);
            //double f(double t) => ((t - 0.25) >= 0 ? 1 : 0) * ((t - 0.75) <= 0 ? 1 : 0);

            double q(double x, double y, double t) => x * y + t + 1;
            double g(double x, double y, double t) => F(x, y, t) - f(t) * q(x, y, t);
            double w(double x, double y, double t) => 1 + x * x + y * y; // w(x, y) only
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

            sfem.Domain(ref node, ref elem, ref dirichlet, xlim, ylim, T, Nx, Ny, Nt);

            // Integral observation
            double[] omega = new double[Nt + 1];
            double[] omega_noise = new double[Nt + 1];
            Random random = new Random();
            // avoid initial condition noise
            for (int nt = 1; nt <= Nt; nt++) {
                omega_noise[nt] = 2 * random.NextDouble() - 1;
            }
            double norm_noise = Math.Sqrt(slvs.NormL2(omega_noise, ht));
            double eps = 0.01;
            for (int nt = 0; nt <= Nt; nt++) {
                double[] temp = new double[(Nx + 1) * (Ny + 1)];
                for (int ny = 0; ny <= Ny; ny++) {
                    for (int nx = 0; nx <= Nx; nx++) {
                        temp[ny * (Nx + 1) + nx] = w(hx * nx, hy * ny, ht * nt) * u(hx * nx, hy * ny, ht * nt);
                    }
                }
                double[] tempx = new double[Nx + 1];
                double[] tempy = new double[Ny + 1];
                for (int ny = 0; ny <= Ny; ny++) {
                    for (int nx = 0; nx <= Nx; nx++) {
                        tempx[nx] = temp[ny * (Nx + 1) + nx];
                    }
                    for (int i = 0; i < Nx; i++) {
                        double funx(double x) => tempx[i] * (1 - x) + tempx[i + 1] * x;
                        tempy[ny] += slvs.Gauss1D(funx) * hx;
                    }
                }
                for (int j = 0; j < Ny; j++) {
                    double funy(double y) => tempy[j] * (1 - y) + tempy[j + 1] * y;
                    omega[nt] += slvs.Gauss1D(funy) * hy;
                }
                omega[nt] = omega[nt] - eps * omega_noise[nt] / norm_noise;
            }
            double[] fh = ifem.CG(node, elem, dirichlet, Nx, Ny, Nt, T, hx, hy, ht, jacobi, omega, gamma, eps, w, axx, ayy, u0, dxu0, dyu0, f, q, g);
            //Console.ReadLine();
        }
    }
}
