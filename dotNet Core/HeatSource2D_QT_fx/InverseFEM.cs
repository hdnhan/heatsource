using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace HeatSource2D_QT_fx {
    class InverseFEM {
        Func<double, double, double> zero2 = (x, y) => 0;
        Func<double, double, double, double> zero3 = (x, y, t) => 0;
        Func<double, double, double> one2 = (x, y) => 1;
        Func<double, double, double, double> one3 = (x, y, t) => 1;

        public double[] CG(double[,] node, int[,] elem, int[] dirichlet, int Nx, int Ny, int Nt, double hx, double hy, double ht, double jacobi, double[] omega, double gamma, double eps,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy, Func<double, double, double> u0, double[] U0,
            Func<double, double, double> f, Func<double, double, double, double> q, Func<double, double, double, double> g) {

            String path_infor = Path.Combine(Environment.CurrentDirectory, @"..\Results\QT_fx_" + Nx + "_" + eps * 100 + "%_log.txt");
            String path_fh = Path.Combine(Environment.CurrentDirectory, @"..\Results\QT_fx_" + Nx + "_" + eps * 100 + "%.txt");
            TextWriter log_infor = new StreamWriter(path_infor);

            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();

            int NoN = node.GetLength(0);

            double[] fh = new double[omega.Length];
            double[] fhold = new double[omega.Length];
            double[] r = new double[omega.Length];
            double[] rold = new double[omega.Length];
            double[] d = new double[omega.Length];

            double errorold, error = 0;
            //fh = oprs.Ones(omega.Length);
            // solution of problem (coz fh = 0)                 
            double[] uh = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, U0, oprs.Ones(omega.Length), g);
            for (int i = 0; i < NoN; i++) {
                uh[i] += u0(node[i, 0], node[i, 1]);
            }

            double[] del_lu = delta_lu(uh, omega, Nx, Ny, Nt);

            for (int iter = 0; iter < 100; iter++) {

                // solution of adjoint problem
                double[] p = AdjointProblem(node, elem, Nx, Ny, Nt, ht, dirichlet, jacobi, axx, ayy, del_lu);

                rold = r;
                r = GradJ(Nx, Ny, Nt, hx, hy, ht, p, q, fh, gamma);

                if (iter > 0) {
                    double beta = slvs.NormL2(r, Nx, Ny, hx, hy) / slvs.NormL2(rold, Nx, Ny, hx, hy);
                    d = oprs.Add(r, oprs.Mul(beta, d));
                } else {
                    d = r;
                }

                double[] udh = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, oprs.Zeros(omega.Length), d, q);
                double[] Ad = delta_lu(udh, oprs.Zeros((Nx + 1) * (Ny + 1)), Nx, Ny, Nt);
                double alpha = slvs.NormL2(r, Nx, Ny, hx, hy) / (slvs.NormL2(Ad, Nx, Ny, hx, hy) + gamma * slvs.NormL2(d, Nx, Ny, hx, hy));

                fhold = fh;
                fh = oprs.Add(fh, oprs.Mul(alpha, d));

                errorold = error;
                error = slvs.ErrorL2(f, fh, Nx, Ny, hx, hy);

                Console.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                Console.WriteLine("J: " + J(Nx, Ny, hx, hy, del_lu, fh, gamma).ToString("e"));
                log_infor.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                log_infor.WriteLine("J: " + J(Nx, Ny, hx, hy, del_lu, fh, gamma).ToString("e"));
                log_infor.WriteLine();

                if (Math.Sqrt(slvs.NormL2(del_lu, Nx, Ny, hx, hy)) < 1.1 * eps) {
                    //break;
                }

                if (iter > 0) {
                    if (error > errorold) {
                        error = errorold;
                        fh = fhold;
                        break;
                    }
                    if (errorold / error < 1 + 1e-3)
                        break;
                }
                del_lu = oprs.Add(del_lu, oprs.Mul(alpha, Ad));
                Console.WriteLine();
            }

            double[] del = new double[(Nx + 1) * (Ny + 1)];
            for (int ny = 0; ny <= Ny; ny++) {
                for (int nx = 0; nx <= Nx; nx++) {
                    del[ny * (Nx + 1) + nx] = fh[ny * (Nx + 1) + nx] - f(nx * hx, ny * hy);
                }
            }

            log_infor.WriteLine("Error in L2: " + slvs.ErrorL2(f, fh, Nx, Ny, hx, hy).ToString("e"));
            log_infor.WriteLine("min and max: " + del.Min().ToString("e") + ", " + del.Max().ToString("e"));
            log_infor.WriteLine("Done!");
            File.WriteAllLines(path_fh, fh.Select(dd => dd.ToString()));
            log_infor.Close();

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


        public double[] AdjointProblem(double[,] node, int[,] elem, int Nx, int Ny, int Nt, double ht, int[] dirichlet, double jacobi,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy, double[] del_lu) {
            Operators oprs = new Operators();
            FEM sfem = new FEM();

            double axxT(double x, double y, double t) => axx(x, y, ht * Nt - t);
            double ayyT(double x, double y, double t) => ayy(x, y, ht * Nt - t);

            double[] p = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axxT, ayyT, del_lu, oprs.Zeros(del_lu.Length), zero3);

            for (int nt = 0; nt <= Nt; nt++) {
                for (int ny = 0; ny <= Ny; ny++) {
                    for (int nx = 0; nx <= Nx; nx++) {
                        p[nt * (Nx + 1) * (Ny + 1) + ny * (Nx + 1) + nx] += del_lu[ny * (Nx + 1) + nx];
                    }
                }
            }
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
