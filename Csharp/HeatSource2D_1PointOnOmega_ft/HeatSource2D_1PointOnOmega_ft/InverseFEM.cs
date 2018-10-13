using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace HeatSource2D_1PointOnOmega_ft {
    class InverseFEM {
        Func<double, double, double> zero2 = (x, y) => 0;
        Func<double, double, double, double> zero3 = (x, y, t) => 0;
        Func<double, double, double> one2 = (x, y) => 1;
        Func<double, double, double, double> one3 = (x, y, t) => 1;

        public double[] CG(double[,] node, int[,] elem, int[] dirichlet, int Nx, int Ny, int Nt, double T,
            double hx, double hy, double ht, double jacobi, double[] omega, double gamma, double eps,
            Func<double, double, double, double> w, Func<double, double, double, double> axx, Func<double, double, double, double> ayy,
            Func<double, double, double> u0, Func<double, double, double> dxu0, Func<double, double, double> dyu0,
            Func<double, double> f, Func<double, double, double, double> q, Func<double, double, double, double> g) {

            String path_infor = Path.Combine(Environment.CurrentDirectory, @"..\..\..\..\Results\1Point_ft_" + Nx + "_" + eps * 100 + "%_log1.txt");
            String path_fh = Path.Combine(Environment.CurrentDirectory, @"..\..\..\..\Results\1Point_ft_" + Nx + "_" + eps * 100 + "%_1.txt");
            TextWriter log_infor = new StreamWriter(path_infor);

            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();

            int NoN = node.GetLength(0);

            double[] fh = new double[Nt + 1];
            double[] fhold = new double[Nt + 1];
            double[] r = new double[Nt + 1];
            double[] rold = new double[Nt + 1];
            double[] d = new double[Nt + 1];

            double errorold, error = 0;

            // solution of problem (coz fh = 0)                 
            double[] uh = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, dxu0, dyu0, oprs.Ones(Nt + 1), g);
            for (int i = 0; i < NoN; i++) {
                uh[i] += u0(node[i, 0], node[i, 1]);
            }

            double[] del_lu = Delta_lu(w, uh, omega, Nx, Ny, Nt, hx, hy, ht);

            for (int iter = 0; iter < 50; iter++) {
                // solution of adjoint problem
                double[] p = AdjointProblem(node, elem, Nx, Ny, Nt, T, dirichlet, jacobi, w, axx, ayy, del_lu);

                rold = r;
                r = GradJ(node, p, q, fh, gamma, Nx, Ny, Ny, hx, hy, ht);
                if (iter > 0) {
                    double beta = slvs.NormL2(r, ht) / slvs.NormL2(rold, ht);
                    d = oprs.Add(r, oprs.Mul(beta, d));
                } else {
                    d = r;
                }

                double[] udh = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, zero2, zero2, d, q);
                double[] Ad = Delta_lu(w, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double alpha = slvs.NormL2(r, ht) / (slvs.NormL2(Ad, ht) + gamma * slvs.NormL2(d, ht));

                fhold = fh;
                fh = oprs.Add(fh, oprs.Mul(alpha, d));

                errorold = error;
                error = slvs.ErrorL2(f, fh, ht);

                Console.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                Console.WriteLine("J: " + J(del_lu, fh, gamma, ht).ToString("e"));
                log_infor.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                log_infor.WriteLine("J: " + J(del_lu, fh, gamma, ht).ToString("e"));
                log_infor.WriteLine();

                if (Math.Sqrt(slvs.NormL2(del_lu, ht)) < 1.1 * eps) {
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

            double[] del = new double[node.GetLength(0)];
            for (int nt = 0; nt <= Nt; nt++) {
                del[nt] = fh[nt] - f(ht * nt);
            }
            log_infor.WriteLine("Error in L2: " + slvs.ErrorL2(f, fh, ht).ToString("e"));
            log_infor.WriteLine("min and max: " + del.Min().ToString("e") + ", " + del.Max().ToString("e"));
            log_infor.WriteLine("Done!");
            File.WriteAllLines(path_fh, fh.Select(dd => dd.ToString()));
            log_infor.Close();

            return fh;
        }


        public double[] GradJ(double[,] node, double[] p, Func<double, double, double, double> q, double[] fh, double gamma,
            int Nx, int Ny, int Nt, double hx, double hy, double ht) {
            Operators oprs = new Operators();
            int NoN = node.GetLength(0);
            double[] res = new double[fh.Length];
            double[] pq = new double[NoN];
            for (int non = 0; non < NoN; non++) {
                pq[non] = p[non] * q(node[non, 0], node[non, 1], node[non, 2]);
            }
            res = Delta_lu(one3, pq, oprs.Zeros(fh.Length), Nx, Ny, Nt, hx, hy, ht);

            return oprs.Add(oprs.Mul(-1, res), oprs.Mul(gamma, fh));
        }


        public double[] AdjointProblem(double[,] node, int[,] elem, int Nx, int Ny, int Nt, double T, int[] dirichlet, double jacobi,
            Func<double, double, double, double> w, Func<double, double, double, double> axx, Func<double, double, double, double> ayy, double[] del_lu) {
            Operators oprs = new Operators();
            FEM sfem = new FEM();

            double axxT(double x, double y, double t) => axx(x, y, T - t);
            double ayyT(double x, double y, double t) => ayy(x, y, T - t);

            double[] p = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axxT, ayyT, zero2, zero2, oprs.flip(del_lu, Nt), w);
            return oprs.flip(p, Nx, Ny, Nt);
        }


        public double J(double[] del_lu, double[] fh, double gamma, double ht) {
            Solvers slvs = new Solvers();
            Operators oprs = new Operators();
            return 0.5 * slvs.NormL2(del_lu, ht) + 0.5 * gamma * slvs.NormL2(fh, ht);
        }


        public double[] Delta_lu(Func<double, double, double, double> w, double[] uh, double[] omega,
            int Nx, int Ny, int Nt, double hx, double hy, double ht) {

            Solvers slvs = new Solvers();
            double[] res = new double[omega.Length];
            for (int nt = 0; nt <= Nt; nt++) {
                double[] temp = new double[(Nx + 1) * (Ny + 1)];
                for (int ny = 0; ny <= Ny; ny++) {
                    for (int nx = 0; nx <= Nx; nx++) {
                        temp[ny * (Nx + 1) + nx] = w(hx * nx, hy * ny, ht * nt) * uh[nt * (Ny + 1) * (Nx + 1) + ny * (Nx + 1) + nx];
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
                    res[nt] += slvs.Gauss1D(funy) * hy;
                }
                res[nt] = res[nt] - omega[nt];
            }
            return res;
        }
    }
}
