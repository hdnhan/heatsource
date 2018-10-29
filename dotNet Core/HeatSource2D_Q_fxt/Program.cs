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

            // Example:
            double u(double x, double y, double t) => Math.Exp(t) * x * (1 - x) * Math.Sin(Math.PI * y);
            double F(double x, double y, double t) => u(x, y, t) + 2 * Math.Exp(t) * Math.Sin(Math.PI * y) + Math.PI * Math.PI * u(x, y, t);

            double u0(double x, double y) => u(x, y, 0);
            double dxu0(double x, double y) => (1 - 2 * x) * Math.Sin(Math.PI * y);
            double dyu0(double x, double y) => Math.PI * x * (1 - x) * Math.Cos(Math.PI * y);

            double f(double x, double y, double t) => Math.Sin(Math.PI * x) * y * (1 - y) * (1 + t * t);
            double q(double x, double y, double t) => x * y + t + 1;
            double g(double x, double y, double t) => F(x, y, t) - f(x, y, t) * q(x, y, t);
            double gamma = 1e-5;

            int nn = 4;
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

            // Observation
            double[] omega = new double[node.GetLength(0)];
            Random random = new Random();
            // avoid dirichlet and initial conditions noise
            for (int nt = 1; nt <= Nt; nt++) {
                for (int ny = 1; ny < Ny; ny++) {
                    for (int nx = 1; nx < Nx; nx++) {
                        omega[nt * (Nx + 1) * (Ny + 1) + ny * (Nx + 1) + nx] = 2 * random.NextDouble() - 1;
                    }
                }
            }

            double norm_noise = Math.Sqrt(slvs.NormL2(omega, elem, jacobi));
            double eps = 0.01;
            for (int i = 0; i < node.GetLength(0); i++) {
                omega[i] = u(node[i, 0], node[i, 1], node[i, 2]) - eps * omega[i] / norm_noise;
            }

            double[] fh = ifem.CG(node, elem, dirichlet, Nx, Ny, Nt, T, jacobi, omega, gamma, eps, axx, ayy, u0, dxu0, dyu0, f, q, g);
        }
    }
}
