using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;


namespace HeatSource2D_Q_fxt {
    class InverseFEM {
        Func<double, double, double> zero2 = (x, y) => 0;
        Func<double, double, double, double> zero3 = (x, y, t) => 0;
        Func<double, double, double> one2 = (x, y) => 1;
        Func<double, double, double, double> one3 = (x, y, t) => 1;

        public double[] CG(double[,] node, int[,] elem, int[] dirichlet, int Nx, int Ny, int Nt, double T, double jacobi, double[] omega, double gamma, double eps,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy,
            Func<double, double, double> u0, Func<double, double, double> dxu0, Func<double, double, double> dyu0,
            Func<double, double, double, double> f, Func<double, double, double, double> q, Func<double, double, double, double> g) {

            String path_infor = Path.Combine(Environment.CurrentDirectory, @"..\..\..\..\Results\Q_fxt_" + Nx + "_" + eps * 100 + "%_log.txt");
            String path_fh = Path.Combine(Environment.CurrentDirectory, @"..\..\..\..\Results\Q_fxt_" + Nx + "_" + eps * 100 + "%.txt");
            TextWriter log_infor = new StreamWriter(path_infor);

            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();

            int NoN = node.GetLength(0);

            double[] fh = new double[NoN];
            double[] fhold = new double[NoN];
            double[] r = new double[NoN];
            double[] rold = new double[NoN];
            double[] d = new double[NoN];

            double errorold, error = 0;

            // solution of problem (coz fh = 0)                 
            double[] uh = sfem.SolveFEM(node, elem, dirichlet, jacobi, axx, ayy, dxu0, dyu0, oprs.Ones(NoN), g);
            for (int i = 0; i < NoN; i++) {
                uh[i] += u0(node[i, 0], node[i, 1]);
            }
            //
            double[] del_lu = oprs.Sub(uh, omega);

            for (int iter = 0; iter < 50; iter++) {
                // solution of adjoint problem
                double[] p = AdjointProblem(node, elem, Nx, Ny, Nt, T, dirichlet, jacobi, axx, ayy, del_lu);

                rold = r;
                r = GradJ(node, p, q, fh, gamma);
                if (iter > 0) {
                    double beta = slvs.NormL2(r, elem, jacobi) / slvs.NormL2(rold, elem, jacobi);
                    d = oprs.Add(r, oprs.Mul(beta, d));
                } else {
                    d = r;
                }

                double[] udh = sfem.SolveFEM(node, elem, dirichlet, jacobi, axx, ayy, zero2, zero2, d, q); // Ad = udh
                double alpha = slvs.NormL2(r, elem, jacobi) / (slvs.NormL2(udh, elem, jacobi) + gamma * slvs.NormL2(d, elem, jacobi));

                fhold = fh;
                fh = oprs.Add(fh, oprs.Mul(alpha, d));

                errorold = error;
                error = slvs.ErrorL2(f, fh, node, elem, jacobi);

                Console.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                Console.WriteLine("J: " + J(elem, jacobi, del_lu, fh, gamma).ToString("e"));

                log_infor.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                log_infor.WriteLine("J: " + J(elem, jacobi, del_lu, fh, gamma).ToString("e"));
                log_infor.WriteLine();

                if (Math.Sqrt(slvs.NormL2(del_lu, elem, jacobi)) < 1.1 * eps) {
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
                del_lu = oprs.Add(del_lu, oprs.Mul(alpha, udh));
                Console.WriteLine();
            }

            double[] del = new double[node.GetLength(0)];
            for (int i = 0; i < node.GetLength(0); i++) {
                del[i] = fh[i] - f(node[i, 0], node[i, 1], node[i, 2]);
            }
            log_infor.WriteLine("Error in L2: " + slvs.ErrorL2(f, fh, node, elem, jacobi).ToString("e"));
            log_infor.WriteLine("min and max: " + del.Min().ToString("e") + ", " + del.Max().ToString("e"));
            log_infor.WriteLine("Done!");
            File.WriteAllLines(path_fh, fh.Select(dd => dd.ToString()));
            log_infor.Close();

            return fh;
        }


        public double[] GradJ(double[,] node, double[] p, Func<double, double, double, double> q, double[] fh, double gamma) {
            int NoN = node.GetLength(0);
            double[] res = new double[NoN];
            for (int non = 0; non < NoN; non++) {
                res[non] = -p[non] * q(node[non, 0], node[non, 1], node[non, 2]) - gamma * fh[non];
            }
            return res;
        }


        public double[] AdjointProblem(double[,] node, int[,] elem, int Nx, int Ny, int Nt, double T, int[] dirichlet, double jacobi,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy, double[] del_lu) {
            Operators oprs = new Operators();
            FEM sfem = new FEM();

            double axxT(double x, double y, double t) => axx(x, y, T - t);
            double ayyT(double x, double y, double t) => ayy(x, y, T - t);

            double[] p = sfem.SolveFEM(node, elem, dirichlet, jacobi, axxT, ayyT, zero2, zero2, oprs.flip(del_lu, Nx, Ny, Nt), one3);
            return oprs.flip(p, Nx, Ny, Nt);
        }


        public double J(int[,] elem, double jacobi, double[] del_lu, double[] fh, double gamma) {
            Solvers slvs = new Solvers();
            Operators oprs = new Operators();
            return 0.5 * slvs.NormL2(del_lu, elem, jacobi) + 0.5 * gamma * slvs.NormL2(fh, elem, jacobi);
        }
    }
}
