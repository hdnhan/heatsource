using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_QT_fx {
    class InverseFEM {
        Func<double, double, double> zero2 = (x, y) => 0;
        Func<double, double, double, double> zero3 = (x, y, t) => 0;
        Func<double, double, double> one2 = (x, y) => 1;
        Func<double, double, double, double> one3 = (x, y, t) => 1;

        public double[] CG(double[,] node, int[,] elem, int[] dirichlet, int Nx, int Ny, int Nt, double hx, double hy, double ht, double jacobi, double[] omega, double gamma, double eps,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy,
            Func<double, double, double> u0, Func<double, double, double> dxu0, Func<double, double, double> dyu0,
            Func<double, double, double> f, Func<double, double, double, double> q, Func<double, double, double, double> g) {

            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();

            int NoN = node.GetLength(0);

            double[] fh = new double[omega.Length];
            double[] fhold = new double[omega.Length];
            double[] r = new double[omega.Length];
            double[] rold = new double[omega.Length];
            double[] d = new double[omega.Length];

            /*
            double[] fe = new double[NoN];
            for (int i = 0; i < NoN; i++) {
                fe[i] = f(node[i, 0], node[i, 1], node[i, 2]);
            }*/

            double errorold, error = 0;

            // solution of problem (coz fh = 0)                 
            double[] uh = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, dxu0, dyu0, oprs.Ones(omega.Length), g);
            for (int i = 0; i < NoN; i++) {
                uh[i] += u0(node[i, 0], node[i, 1]);
            }
            //uh = oprs.Add(uh, sfem.SolveFEM(node, elem, dirichlet, jacobi, axx, ayy, zero2, zero2, fe, q));
            double[] del_lu = delta_lu(uh, omega, Nx, Ny, Nt);
            //Console.WriteLine("del_lu min and max: " + del_lu.Min() + ", " + del_lu.Max());

            for (int iter = 0; iter < 50; iter++) {
                // solution of adjoint problem
                double[] p = AdjointProblem(node, elem, Nx, Ny, Nt, dirichlet, jacobi, axx, ayy, del_lu);

                rold = r;
                r = GradJ(Nx, Ny, Nt, hx, hy, ht, p, q, fh, gamma);
                if (iter > 0) {
                    double beta = slvs.NormL2(r, Nx, Ny, hx, hy) / slvs.NormL2(rold, Nx, Ny, hx, hy);
                    d = oprs.Add(r, oprs.Mul(beta, d));
                } else {
                    d = r;
                }

                double[] udh = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, zero2, zero2, d, q);
                double[] Ad = delta_lu(udh, oprs.Zeros((Nx + 1) * (Ny + 1)), Nx, Ny, Nt);
                double alpha = slvs.NormL2(r, Nx, Ny, hx, hy) / (slvs.NormL2(Ad, Nx, Ny, hx, hy) + gamma * slvs.NormL2(d, Nx, Ny, hx, hy));

                fhold = fh;
                fh = oprs.Add(fh, oprs.Mul(alpha, d));

                double[] del = new double[omega.GetLength(0)];
                for (int ny = 0; ny <= Ny; ny++) {
                    for (int nx = 0; nx <= Nx; nx++) {
                        del[ny * (Nx + 1) + nx] = fh[ny * (Nx + 1) + nx] - f(nx * hx, ny * hy);
                        //del[ny * (Nx + 1) + nx] = test1(nx * hx, ny * hy);
                    }
                }
                errorold = error;
                error = Math.Sqrt(slvs.NormL2(del, Nx, Ny, hx, hy));

                Console.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                Console.WriteLine("J: " + J(Nx, Ny, hx, hy, del_lu, fh, gamma).ToString("e"));

                if (Math.Sqrt(slvs.NormL2(del_lu, Nx, Ny, hx, hy)) < 1.1 * eps) {
                    //break;
                }

                if (iter > 0) {
                    if (error > errorold) {
                        error = errorold;
                        fh = fhold;
                        break;
                    }
                }
                del_lu = oprs.Add(del_lu, oprs.Mul(alpha, udh));
                Console.WriteLine();
            }
            return fh;
        }


        public double[] GradJ(int Nx, int Ny, int Nt, double hx, double hy, double ht, 
            double[] p, Func<double, double, double, double> q, double[] fh, double gamma) {
            Solvers slvs = new Solvers();

            double[] res = new double[(Nx + 1) * (Ny + 1)];
            double[] temp = new double[Nt + 1];
            for (int ny = 0; ny <= Ny; ny++) {
                for (int nx = 0; nx <= Nx; nx++) {
                    for (int nt = 0; nt <= Nt; nt++) {
                        temp[nt] = p[nt * (Nx + 1) * (Ny + 1) + ny * (Nx + 1) + nx] * q(nx * hx, ny * hy, nt * ht);
                    }
                    for (int i = 0; i < Nt; i++) {
                        double fun(double t) => temp[i] * (1 - t) + temp[i + 1] * t;
                        res[ny * (Nx + 1) + nx] += slvs.Gauss1D(fun) * ht;
                    }
                    res[ny * (Nx + 1) + nx] = -res[ny * (Nx + 1) + nx] - gamma * fh[ny * (Nx + 1) + nx];
                }
            }
            return res;
        }


        public double[] AdjointProblem(double[,] node, int[,] elem, int Nx, int Ny, int Nt, int[] dirichlet, double jacobi,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy, double[] del_lu) {
            Operators oprs = new Operators();
            FEM sfem = new FEM();
            double[] p = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, zero2, zero2, oprs.flip(del_lu, Nt), one3);
            return oprs.flip(p, Nx, Ny, Nt);
        }


        public double J(int Nx, int Ny, double hx, double hy, double[] del_lu, double[] fh, double gamma) {
            Solvers slvs = new Solvers();
            Operators oprs = new Operators();
            return 0.5 * slvs.NormL2(del_lu, Nx, Ny, hx, hy) + 0.5 * gamma * slvs.NormL2(fh, Nx, Ny, hx, hy);
        }

        public double[] delta_lu(double[] uh, double[] omega, int Nx, int Ny, int Nt) {
            double[] res = new double[omega.Length];
            for (int i = 0; i < omega.Length; i++) {
                res[i] = uh[Nt * (Nx + 1) * (Ny + 1) + i] - omega[i];
            }
            return res;
        }
    }
}
