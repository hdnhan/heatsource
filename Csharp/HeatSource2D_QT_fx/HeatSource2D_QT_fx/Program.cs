using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_QT_fx {
    class Program {
        static void Main(string[] args) {
            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();
            InverseFEM ifem = new InverseFEM();

            // DATA
            double axx(double x, double y, double t) => 1;
            double ayy(double x, double y, double t) => 1;

            // Example:
            double u(double x, double y, double t) => Math.Exp(t) * x * (1 - x) * Math.Sin(Math.PI * y);
            double F(double x, double y, double t) => u(x, y, t) + 2 * Math.Exp(t) * Math.Sin(Math.PI * y) + Math.PI * Math.PI * u(x, y, t);

            double u0(double x, double y) => u(x, y, 0);
            //double dxu0(double x, double y) => (1 - 2 * x) * Math.Sin(Math.PI * y);
            //double dyu0(double x, double y) => Math.PI * x * (1 - x) * Math.Cos(Math.PI * y);

            double f(double x, double y) => Math.Sin(Math.PI * x) * y * (1 - y);
            double q(double x, double y, double t) => x * y + t + 1;
            double g(double x, double y, double t) => F(x, y, t) - f(x, y) * q(x, y, t);
            double gamma = 1e-5;

            int nn = 64;
            int Nx = nn, Ny = nn, Nt = nn;
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

            sfem.Domain(ref node, ref elem, ref dirichlet, xlim, ylim, T, Nx, Ny, Nt);

            // Final bservation (eps => noise) and discrete initial condition
            double[] omega = new double[(Nx + 1) * (Ny + 1)];
            double[] U0 = new double[(Nx + 1) * (Ny + 1)];
            Random random = new Random();
            // avoid dirichlet condition noise
            for (int ny = 1; ny < Ny; ny++) {
                for (int nx = 1; nx < Nx; nx++) {
                    omega[ny * (Nx + 1) + nx] = 2 * random.NextDouble() - 1; // eps * [-1, 1]
                }
            }
            double norm_noise = Math.Sqrt(slvs.NormL2(omega, Nx, Ny, hx, hy));
            double eps = 0.01;
            for (int ny = 0; ny <= Ny; ny++) {
                for (int nx = 0; nx <= Nx; nx++) {
                    omega[ny * (Nx + 1) + nx] = u(nx * hx, ny * hy, T) - eps * omega[ny * (Nx + 1) + nx] / norm_noise;
                    U0[ny * (Nx + 1) + nx] = u0(nx * hx, ny * hy);
                }
            }
            
            double[] fh = ifem.CG(node, elem, dirichlet, Nx, Ny, Nt, hx, hy, ht, jacobi, omega, gamma, eps, axx, ayy, u0, U0, f, q, g);
            //Console.ReadLine();
        }
    }
}
