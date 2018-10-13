using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_Q_fxt {
    class Solvers {

        // Tính tích phân 3 lớp f(x, y) x: 0 -> 1, y: 0 -> 1-x, z: 1-x-y
        public double Gauss3D(Func<double, double, double, double> f3) {
            double[] w = new double[] { 0.277777777777778, 0.444444444444444, 0.277777777777778 };
            double[] p = new double[] { 0.112701665379258, 0.5, 0.887298334620742 };

            double f2(double x, double y) => w[0] * (1 - x - y) * f3(x, y, (1 - x - y) * p[0]) +
                                             w[1] * (1 - x - y) * f3(x, y, (1 - x - y) * p[1]) +
                                             w[2] * (1 - x - y) * f3(x, y, (1 - x - y) * p[2]);

            double f1(double x) => w[0] * (1 - x) * f2(x, (1 - x) * p[0]) +
                                   w[1] * (1 - x) * f2(x, (1 - x) * p[1]) +
                                   w[2] * (1 - x) * f2(x, (1 - x) * p[2]);

            return w[0] * f1(p[0]) + w[1] * f1(p[1]) + w[2] * f1(p[2]);
        }

        // Tính tích phân 2 lớp f(x, y) x: 0 -> 1, y: 0 -> 1-x
        public double Gauss2D(Func<double, double, double> f2) {
            double[] w = new double[] { 0.277777777777778, 0.444444444444444, 0.277777777777778 };
            double[] p = new double[] { 0.112701665379258, 0.5, 0.887298334620742 };

            double f1(double x) => w[0] * (1 - x) * f2(x, (1 - x) * p[0]) +
                                   w[1] * (1 - x) * f2(x, (1 - x) * p[1]) +
                                   w[2] * (1 - x) * f2(x, (1 - x) * p[2]);

            return w[0] * f1(p[0]) + w[1] * f1(p[1]) + w[2] * f1(p[2]);
        }

        public double Gauss1D(Func<double, double> f1) {
            double[] w = new double[] { 0.277777777777778, 0.444444444444444, 0.277777777777778 };
            double[] p = new double[] { 0.112701665379258, 0.5, 0.887298334620742 };

            return w[0] * f1(p[0]) + w[1] * f1(p[1]) + w[2] * f1(p[2]);
        }

        public double[] BiCGM(Dictionary<Tuple<int, int>, double> A, double[] b) {

            Operators ops = new Operators();
            int leng = b.Length;
            double[] x = new double[leng];
            double[] r1 = new double[leng];
            double[] r2 = new double[leng];
            double[] p1 = new double[leng];
            double[] p2 = new double[leng];
            double[] q1 = new double[leng];
            double[] q2 = new double[leng];
            double res, resold = 0, beta, alpha;
            double resfix = ops.InnerProduct(b, b);

            r1 = b; //DOs.Sub(b, DOs.Mul(A, x)); // r1 = b-A*x0
            r2 = r1;
            for (int iter = 0; iter <= leng; iter++) {
                res = ops.InnerProduct(r1, r2); // res = r1'*r2

                if (Math.Abs(res) <= 1e-12 * resfix) {
                    break;
                }
                if (iter > 0) {
                    beta = res / resold;
                    p1 = ops.Add(r1, ops.Mul(beta, p1)); // p1 = r1+ beta*p1
                    p2 = ops.Add(r2, ops.Mul(beta, p2));
                } else {
                    p1 = r1;
                    p2 = r2;
                }

                q1 = ops.Mul(A, p1); // q1 = A * p1
                q2 = ops.MulI(A, p2); // q2=A'*p2
                alpha = res / ops.InnerProduct(p2, q1); //alpha = res / ( p2' * q1 )
                x = ops.Add(x, ops.Mul(alpha, p1)); // x = x + alpha * p1
                r1 = ops.Sub(r1, ops.Mul(alpha, q1)); // r1 = r1 - alpha * q1
                r2 = ops.Sub(r2, ops.Mul(alpha, q2));
                resold = res;
            }
            return x;
        }

        // All error function in L2
        // f(x, y, t)
        public double ErrorL2(Func<double, double, double, double> f, double[] fh, double[,] node, int[,] elem, double jacobi) {
            double res = 0;
            int NoE = elem.GetLength(0);
            for (int noe = 0; noe < NoE; noe++) {
                double[] x = new double[] { node[elem[noe, 0], 0], node[elem[noe, 1], 0], node[elem[noe, 2], 0], node[elem[noe, 3], 0] };
                double[] y = new double[] { node[elem[noe, 0], 1], node[elem[noe, 1], 1], node[elem[noe, 2], 1], node[elem[noe, 3], 1] };
                double[] t = new double[] { node[elem[noe, 0], 2], node[elem[noe, 1], 2], node[elem[noe, 2], 2], node[elem[noe, 3], 2] };

                double[] del = new double[] { fh[elem[noe, 0]] - f(x[0], y[0], t[0]), fh[elem[noe, 1]] - f(x[1], y[1], t[1]),
                                              fh[elem[noe, 2]] - f(x[2], y[2], t[2]), fh[elem[noe, 3]] - f(x[3], y[3], t[3]) };

                double del_hat(double xi, double eta, double zeta) => del[0] * (1 - xi - eta - zeta) + del[1] * xi + del[2] * eta + del[3] * zeta;
                double del_hat2(double xi, double eta, double zeta) => del_hat(xi, eta, zeta) * del_hat(xi, eta, zeta);
                res += Gauss3D(del_hat2) * jacobi;
            }
            return Math.Sqrt(res);
        }

        // f(x, y)
        public double ErrorL2(Func<double, double, double> f, double[] fh, int Nx, int Ny, double hx, double hy) {
            double res = 0;
            double[] delx2 = new double[Nx + 1];
            double[] dely = new double[Ny + 1];

            for (int ny = 0; ny <= Ny; ny++) {
                for (int nx = 0; nx <= Nx; nx++) {
                    delx2[nx] = (fh[ny * (Nx + 1) + nx] - f(nx * hx, ny * hy)) * (fh[ny * (Nx + 1) + nx] - f(nx * hx, ny * hy));
                }
                for (int i = 0; i < Nx; i++) {
                    double funx(double x) => delx2[i] * (1 - x) + delx2[i + 1] * x;
                    dely[ny] += Gauss1D(funx) * hx;
                }
            }
            for (int j = 0; j < Ny; j++) {
                double funy(double y) => dely[j] * (1 - y) + dely[j + 1] * y;
                res += Gauss1D(funy) * hy;
            }
            return Math.Sqrt(res);
        }

        // f(t)
        public double ErrorL2(Func<double, double> f, double[] fh, double ht) {
            double res = 0;
            int Nt = fh.Length - 1;
            double[] del2 = new double[Nt + 1];

            for (int nt = 0; nt <= Nt; nt++) {
                del2[nt] = (fh[nt] - f(nt * ht)) * (fh[nt] - f(nt * ht));
            }

            for (int nt = 0; nt < Nt; nt++) {
                double fun(double t) => del2[nt] * (1 - t) + del2[nt + 1] * t;
                res += Gauss1D(fun) * ht;
            }

            return Math.Sqrt(res);
        }


        // All Norm in L2
        // f(x, y, t)
        public double NormL2(double[] fh, int[,] elem, double jacobi) {
            double res = 0;
            int NoE = elem.GetLength(0);
            for (int noe = 0; noe < NoE; noe++) {
                double[] z = new double[] { fh[elem[noe, 0]], fh[elem[noe, 1]], fh[elem[noe, 2]], fh[elem[noe, 3]] };
                double fh_hat2(double xi, double eta, double zeta)
                    => z[0] * z[0] * (1 - xi - eta - zeta) + z[1] * z[1] * xi + z[2] * z[2] * eta + z[3] * z[3] * zeta;
                res += Gauss3D(fh_hat2) * jacobi;

            }
            return res;
        }

        // f(x, y)
        public double NormL2(double[] uh, int Nx, int Ny, double hx, double hy) {
            double res = 0;
            double[] tempx = new double[Nx + 1];
            double[] tempy = new double[Ny + 1];

            for (int ny = 0; ny <= Ny; ny++) {
                for (int nx = 0; nx <= Nx; nx++) {
                    tempx[nx] = uh[ny * (Nx + 1) + nx] * uh[ny * (Nx + 1) + nx];
                }
                for (int i = 0; i < Nx; i++) {
                    double funx(double x) => tempx[i] * (1 - x) + tempx[i + 1] * x;
                    tempy[ny] += Gauss1D(funx) * hx;
                }
            }
            for (int j = 0; j < Ny; j++) {
                double funy(double y) => tempy[j] * (1 - y) + tempy[j + 1] * y;
                res += Gauss1D(funy) * hy;
            }
            return res;
        }

        // f(t)
        public double NormL2(double[] fh, double ht) {
            double res = 0;
            int Nt = fh.Length - 1;
            for (int nt = 0; nt < Nt; nt++) {
                double fun(double t) => fh[nt] * fh[nt] * (1 - t) + fh[nt + 1] * fh[nt + 1] * t;
                res += Gauss1D(fun) * ht;
            }
            return res;
        }
    }
}
