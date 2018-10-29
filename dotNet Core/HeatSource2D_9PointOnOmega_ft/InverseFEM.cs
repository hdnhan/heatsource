using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace HeatSource2D_9PointOnOmega_ft {
    class InverseFEM {
        Func<double, double, double> zero2 = (x, y) => 0;
        Func<double, double, double, double> zero3 = (x, y, t) => 0;
        Func<double, double, double> one2 = (x, y) => 1;
        Func<double, double, double, double> one3 = (x, y, t) => 1;

        public double[] CG(double[,] node, int[,] elem, int[] dirichlet, int Nx, int Ny, int Nt, double T,
            double hx, double hy, double ht, double jacobi, double gamma, double eps,
            double[] omega0, double[] omega1, double[] omega2, double[] omega3, double[] omega4, double[] omega5, double[] omega6, double[] omega7, double[] omega8,
            Func<double, double, double, double> w0, Func<double, double, double, double> w1, Func<double, double, double, double> w2,
            Func<double, double, double, double> w3, Func<double, double, double, double> w4, Func<double, double, double, double> w5,
            Func<double, double, double, double> w6, Func<double, double, double, double> w7, Func<double, double, double, double> w8, 
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy,
            Func<double, double, double> u0, Func<double, double, double> dxu0, Func<double, double, double> dyu0,
            Func<double, double, double> f, Func<double, double, double, double> q, Func<double, double, double, double> g) {

            String path_infor = Path.Combine(Environment.CurrentDirectory, @"..\Results\9Point_ft_" + Nx + "_" + eps * 100 + "%_log.txt");
            String path_fh = Path.Combine(Environment.CurrentDirectory, @"..\Results\9Point_ft_" + Nx + "_" + eps * 100 + "%.txt");
            TextWriter log_infor = new StreamWriter(path_infor);

            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();

            int NoN = node.GetLength(0);

            double[] fh = new double[(Nx+1)*(Ny+1)];
            double[] fhold = new double[(Nx + 1) * (Ny + 1)];
            double[] r = new double[(Nx + 1) * (Ny + 1)];
            double[] rold = new double[(Nx + 1) * (Ny + 1)];
            double[] d = new double[(Nx + 1) * (Ny + 1)];

            double errorold, error = 0;

            // solution of problem (coz fh = 0)                 
            double[] uh = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, dxu0, dyu0, oprs.Ones(Nt + 1), g,
                oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3,
                oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3);
            for (int i = 0; i < NoN; i++) {
                uh[i] += u0(node[i, 0], node[i, 1]);
            }

            double[] del_lu0 = Delta_lu(w0, uh, omega0, Nx, Ny, Nt, hx, hy, ht);
            double[] del_lu1 = Delta_lu(w1, uh, omega1, Nx, Ny, Nt, hx, hy, ht);
            double[] del_lu2 = Delta_lu(w2, uh, omega2, Nx, Ny, Nt, hx, hy, ht);
            double[] del_lu3 = Delta_lu(w3, uh, omega3, Nx, Ny, Nt, hx, hy, ht);
            double[] del_lu4 = Delta_lu(w4, uh, omega4, Nx, Ny, Nt, hx, hy, ht);
            double[] del_lu5 = Delta_lu(w5, uh, omega5, Nx, Ny, Nt, hx, hy, ht);
            double[] del_lu6 = Delta_lu(w6, uh, omega6, Nx, Ny, Nt, hx, hy, ht);
            double[] del_lu7 = Delta_lu(w7, uh, omega7, Nx, Ny, Nt, hx, hy, ht);
            double[] del_lu8 = Delta_lu(w8, uh, omega8, Nx, Ny, Nt, hx, hy, ht);

            for (int iter = 0; iter < 50; iter++) {
                // solution of adjoint problem
                double[] p = AdjointProblem(node, elem, Nx, Ny, Nt, T, dirichlet, jacobi, w0, w1, w2, w3, w4, w5, w6, w7, w8, 
                    axx, ayy, del_lu0, del_lu1, del_lu2, del_lu3, del_lu4, del_lu5, del_lu6, del_lu7, del_lu8);

                rold = r;
                r = GradJ(node, p, q, fh, gamma, Nx, Ny, Ny, hx, hy, ht);
                if (iter > 0) {
                    double beta = slvs.NormL2(r, ht) / slvs.NormL2(rold, ht);
                    d = oprs.Add(r, oprs.Mul(beta, d));
                } else {
                    d = r;
                }

                double[] udh = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axx, ayy, zero2, zero2, d, q,
                    oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3,
                    oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3, oprs.Zeros(Nt + 1), zero3);

                double[] Ad0 = Delta_lu(w0, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double[] Ad1 = Delta_lu(w1, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double[] Ad2 = Delta_lu(w2, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double[] Ad3 = Delta_lu(w3, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double[] Ad4 = Delta_lu(w4, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double[] Ad5 = Delta_lu(w5, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double[] Ad6 = Delta_lu(w6, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double[] Ad7 = Delta_lu(w7, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);
                double[] Ad8 = Delta_lu(w8, udh, oprs.Zeros(Nt + 1), Nx, Ny, Nt, hx, hy, ht);

                double alpha = slvs.NormL2(r, ht) / (slvs.NormL2(Ad0, ht) + slvs.NormL2(Ad2, ht) + slvs.NormL2(Ad3, ht) +
                                                     slvs.NormL2(Ad3, ht) + slvs.NormL2(Ad4, ht) + slvs.NormL2(Ad5, ht) +
                                                     slvs.NormL2(Ad6, ht) + slvs.NormL2(Ad7, ht) + slvs.NormL2(Ad8, ht) +
                                                     gamma * slvs.NormL2(d, ht));

                fhold = fh;
                fh = oprs.Add(fh, oprs.Mul(alpha, d));

                errorold = error;
                error = slvs.ErrorL2(f, fh, Nx, Ny, hx, hy);

                Console.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                Console.WriteLine("J: " + J(del_lu0, del_lu1, del_lu2, del_lu3, del_lu4, del_lu5, del_lu6, del_lu7, del_lu8, fh, gamma, ht).ToString("e"));
                log_infor.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                log_infor.WriteLine("J: " + J(del_lu0, del_lu1, del_lu2, del_lu3, del_lu4, del_lu5, del_lu6, del_lu7, del_lu8, fh, gamma, ht).ToString("e"));
                log_infor.WriteLine();

                //if (Math.Sqrt(slvs.NormL2(del_lu, ht)) < 1.1 * eps) {
                    //break;
                //}

                if (iter > 0) {
                    if (error > errorold) {
                        error = errorold;
                        fh = fhold;
                        break;
                    }
                    if (errorold / error < 1 + 1e-3)
                        break;
                }
                del_lu0 = oprs.Add(del_lu0, oprs.Mul(alpha, Ad0));
                del_lu1 = oprs.Add(del_lu1, oprs.Mul(alpha, Ad1));
                del_lu2 = oprs.Add(del_lu2, oprs.Mul(alpha, Ad2));
                del_lu3 = oprs.Add(del_lu3, oprs.Mul(alpha, Ad3));
                del_lu4 = oprs.Add(del_lu4, oprs.Mul(alpha, Ad4));
                del_lu5 = oprs.Add(del_lu5, oprs.Mul(alpha, Ad5));
                del_lu6 = oprs.Add(del_lu6, oprs.Mul(alpha, Ad6));
                del_lu7 = oprs.Add(del_lu7, oprs.Mul(alpha, Ad7));
                del_lu8 = oprs.Add(del_lu8, oprs.Mul(alpha, Ad8));
                Console.WriteLine();
            }

            double[] del = new double[(Nx + 1) * (Ny + 1)];
            for (int ny = 0; ny <= Ny; ny++) {
                for (int nx = 0; nx <= Nx; nx++) {
                    del[ny * (Ny + 1) + nx] = fh[ny * (Ny + 1) + nx] - f(hx * nx, hy * ny);
                }
            }
            log_infor.WriteLine("Error in L2: " + slvs.ErrorL2(f, fh, Nx, Ny, hx, hy).ToString("e"));
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
            Func<double, double, double, double> w0, Func<double, double, double, double> w1, Func<double, double, double, double> w2,
            Func<double, double, double, double> w3, Func<double, double, double, double> w4, Func<double, double, double, double> w5,
            Func<double, double, double, double> w6, Func<double, double, double, double> w7, Func<double, double, double, double> w8,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy, double[] del_lu0, double[] del_lu1, 
            double[] del_lu2, double[] del_lu3, double[] del_lu4, double[] del_lu5, double[] del_lu6, double[] del_lu7, double[] del_lu8) {
            Operators oprs = new Operators();
            FEM sfem = new FEM();

            double axxT(double x, double y, double t) => axx(x, y, T - t);
            double ayyT(double x, double y, double t) => ayy(x, y, T - t);

            double[] p = sfem.SolveFEM(node, elem, dirichlet, Nx, Ny, jacobi, axxT, ayyT, zero2, zero2, 
                oprs.flip(del_lu0, Nt), w0, oprs.flip(del_lu2, Nt), w2, oprs.flip(del_lu3, Nt), w3,
                oprs.flip(del_lu3, Nt), w3, oprs.flip(del_lu4, Nt), w4, oprs.flip(del_lu5, Nt), w5,
                oprs.flip(del_lu6, Nt), w6, oprs.flip(del_lu7, Nt), w7, oprs.flip(del_lu8, Nt), w8);
            return oprs.flip(p, Nx, Ny, Nt);
        }


        public double J(double[] del_lu0, double[] del_lu1, double[] del_lu2,
                        double[] del_lu3, double[] del_lu4, double[] del_lu5,
                        double[] del_lu6, double[] del_lu7, double[] del_lu8,
                        double[] fh, double gamma, double ht) {
            Solvers slvs = new Solvers();
            Operators oprs = new Operators();
            return 0.5 * slvs.NormL2(del_lu0, ht) + 0.5 * slvs.NormL2(del_lu1, ht) + 0.5 * slvs.NormL2(del_lu2, ht) +
                   0.5 * slvs.NormL2(del_lu3, ht) + 0.5 * slvs.NormL2(del_lu4, ht) + 0.5 * slvs.NormL2(del_lu5, ht) +
                   0.5 * slvs.NormL2(del_lu6, ht) + 0.5 * slvs.NormL2(del_lu7, ht) + 0.5 * slvs.NormL2(del_lu8, ht) +
                   0.5 * gamma * slvs.NormL2(fh, ht);
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
