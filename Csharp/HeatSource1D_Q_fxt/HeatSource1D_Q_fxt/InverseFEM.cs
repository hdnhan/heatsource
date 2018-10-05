using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource1D_Q_fxt {
    class InverseFEM {
        Func<double, double> zero1 = (x) => 0;
        Func<double, double, double> zero2 = (x, t) => 0;
        Func<double, double> one1 = (x) => 1;
        Func<double, double, double> one2 = (x, t) => 1;

        public double[] CG(double[,] node, int[,] elem, int[] dirichlet, int Nx, int Nt, double jacobi,
            double[] omega, double gamma, double eps, Func<double, double, double> a, Func<double, double> u0, Func<double, double> dxu0,
            Func<double, double, double> f, Func<double, double, double> q, Func<double, double, double> g) {

            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();

            int NoN = node.GetLength(0);

            double[] fh = new double[NoN];
            double[] fhold = new double[NoN];
            double[] r = new double[NoN];
            double[] rold = new double[NoN];
            double[] d = new double[NoN];

            double[] fe = new double[NoN];
            for (int i = 0; i < NoN; i++) {
                fe[i] = f(node[i, 0], node[i, 1]);
            }

            double errorold, error = 0;

            // solution of problem (coz fh = 0)           
            double[] uh = sfem.SolveFEM(node, elem, dirichlet, jacobi, a, dxu0, oprs.Ones(NoN), g);
            //fh = fe;
            //uh = oprs.Add(uh, sfem.SolveFEM(node, elem, dirichlet, jacobi, a, zero1, fh, q));
            for (int i = 0; i < NoN; i++) {
                uh[i] += u0(node[i, 0]);
            }
            double[] del_lu = oprs.Sub(uh, omega);
            //oprs.Print(del_lu);
            //oprs.Print(oprs.flip(del_lu, Nx, Nt));
            //Console.WriteLine("del_lu min and max: " + del_lu.Min() + ", " + del_lu.Max());


            for (int iter = 0; iter < 50; iter++) {
                // solution of adjoint problem
                //oprs.Print(del_lu);
                double[] p = AdjointProblem(node, elem, Nx, Nt, dirichlet, jacobi, a, del_lu);
                //Console.WriteLine("p min and max: " + p.Min() + ", " + p.Max());
                //oprs.Print(del_lu);
                // r = - Gradient J(f)
                rold = r;
                r = GradJ(node, p, q, fh, gamma);
                if (iter > 0) {
                    double beta = slvs.NormL2(r, elem, jacobi) / slvs.NormL2(rold, elem, jacobi);
                    d = oprs.Add(r, oprs.Mul(beta, d));
                } else {
                    d = r;
                }

                double[] udh = sfem.SolveFEM(node, elem, dirichlet, jacobi, a, zero1, d, q); // Ad = udh
                double alpha = slvs.NormL2(r, elem, jacobi) / (slvs.NormL2(udh, elem, jacobi) + gamma * slvs.NormL2(d, elem, jacobi));
                //Console.WriteLine("alpha: " + alpha);
                //Console.WriteLine(slvs.NormL2(r, elem, jacobi));
                //Console.WriteLine(slvs.NormL2(udh, elem, jacobi));
                //Console.WriteLine(gamma * slvs.NormL2(d, elem, jacobi));

                fhold = fh;
                fh = oprs.Add(fh, oprs.Mul(alpha, d));

                errorold = error;
                error = slvs.ErrorL2(f, fh, node, elem, jacobi);

                Console.WriteLine(iter + ": ErrorL2: " + error.ToString("e") + ", DeltaError: " + (errorold - error).ToString("e"));
                Console.WriteLine("J: " + J(elem, jacobi, del_lu, fh, gamma).ToString("e"));

                if (Math.Sqrt(slvs.NormL2(del_lu, elem, jacobi)) < 1.1 * eps) {
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
                //Console.WriteLine("del_lu min and max: " + del_lu.Min() + ", " + del_lu.Max());
                //double[] del_f = oprs.Sub(fh, fe);
                //Console.WriteLine("del_f min and max: " + del_f.Min() + ", " + del_f.Max());
                Console.WriteLine();
            }
            return fh;
        }


        public double[] GradJ(double[,] node, double[] p, Func<double, double, double> q, double[] fh, double gamma) {
            int NoN = node.GetLength(0);
            double[] res = new double[NoN];
            for (int non = 0; non < NoN; non++) {
                res[non] = -p[non] * q(node[non, 0], node[non, 1]) - gamma * fh[non];
            }

            return res;
        }


        public double[] AdjointProblem(double[,] node, int[,] elem, int Nx, int Nt,
            int[] dirichlet, double jacobi, Func<double, double, double> a, double[] del_lu) {
            Operators oprs = new Operators();
            FEM sfem = new FEM();
            double[] p = sfem.SolveFEM(node, elem, dirichlet, jacobi, a, zero1, oprs.flip(del_lu, Nx, Nt), one2);
            //Console.WriteLine("p min and max: " + p.Min() + ", " + p.Max());
            return oprs.flip(p, Nx, Nt);
        }


        public double J(int[,] elem, double jacobi, double[] del_lu, double[] fh, double gamma) {
            Solvers slvs = new Solvers();
            Operators oprs = new Operators();
            return 0.5 * slvs.NormL2(del_lu, elem, jacobi) + 0.5 * gamma * slvs.NormL2(fh, elem, jacobi);
        }
    }
}
