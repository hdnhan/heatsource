using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpaceTime2D
{
    class Solvers
    {
        // Tính tích phân 3 lớp f(x, y) x: 0 -> 1, y: 0 -> 1-x, z: 1-x-y
        public double Gauss3D(Func<double, double, double, double> f)
        {
            double[] w = new double[] { 0.277777777777778, 0.444444444444444, 0.277777777777778 };
            double[] p = new double[] { 0.112701665379258, 0.5, 0.887298334620742 };

            Func<double, double, double> f2 = (x, y) => w[0] * (1 - x - y) * f(x, y, (1 - x - y) * p[0]) +
                                                        w[1] * (1 - x - y) * f(x, y, (1 - x - y) * p[1]) +
                                                        w[2] * (1 - x - y) * f(x, y, (1 - x - y) * p[2]);

            Func<double, double> f1 = (x) => w[0] * (1 - x) * f2(x, (1 - x) * p[0]) +
                                             w[1] * (1 - x) * f2(x, (1 - x) * p[1]) +
                                             w[2] * (1 - x) * f2(x, (1 - x) * p[2]);

            return w[0] * f1(p[0]) + w[1] * f1(p[1]) + w[2] * f1(p[2]);
        }

        // Tính tích phân 2 lớp f(x, y) x: 0 -> 1, y: 0 -> 1-x
        public double Gauss2D(Func<double, double, double> f)
        {
            double[] w = new double[] { 0.277777777777778, 0.444444444444444, 0.277777777777778 };
            double[] p = new double[] { 0.112701665379258, 0.5, 0.887298334620742 };

            //Func<double, double, double> ff = (x, y) => (1 - x) * f(x, (1 - x) * y);
            //Func<double, double> fff = (x) => w[0] * ff(x, p[0]) + w[1] * ff(x, p[1]) + w[2] * ff(x, p[2]);

            Func<double, double> g = (x) => w[0] * (1 - x) * f(x, (1 - x) * p[0]) +
                                            w[1] * (1 - x) * f(x, (1 - x) * p[1]) +
                                            w[2] * (1 - x) * f(x, (1 - x) * p[2]);
            return w[0] * g(p[0]) + w[1] * g(p[1]) + w[2] * g(p[2]);
        }

        public double Gauss1D(Func<double, double> f)
        {
            double[] w = new double[] { 0.277777777777778, 0.444444444444444, 0.277777777777778 };
            double[] p = new double[] { 0.112701665379258, 0.5, 0.887298334620742 };
            return w[0] * f(p[0]) + w[1] * f(p[1]) + w[2] * f(p[2]);
        }

        public double[] BiCGM(Dictionary<Tuple<int, int>, double> A, double[] b)
        {

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
            for (int iter = 0; iter <= leng; iter++)
            {
                res = ops.InnerProduct(r1, r2); // res = r1'*r2

                if (Math.Abs(res) <= 1e-12 * resfix)
                {
                    break;
                }
                if (iter > 0)
                {
                    beta = res / resold;
                    p1 = ops.Add(r1, ops.Mul(beta, p1)); // p1 = r1+ beta*p1
                    p2 = ops.Add(r2, ops.Mul(beta, p2));
                }
                else
                {
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


        public double ErrorL2(Func<double, double, double, double> u, double[] uh, double[,] node, int[,] elem, double jacobi)
        {
            double res = 0;
            int NoE = elem.GetLength(0);
            for (int noe = 0; noe < NoE; noe++)
            {
                double[] x = new double[] { node[elem[noe, 0], 0], node[elem[noe, 1], 0], node[elem[noe, 2], 0], node[elem[noe, 3], 0] };
                double[] y = new double[] { node[elem[noe, 0], 1], node[elem[noe, 1], 1], node[elem[noe, 2], 1], node[elem[noe, 3], 1] };
                double[] t = new double[] { node[elem[noe, 0], 2], node[elem[noe, 1], 2], node[elem[noe, 2], 2], node[elem[noe, 3], 2] };

                double[] z = new double[] { uh[elem[noe, 0]], uh[elem[noe, 1]], uh[elem[noe, 2]], uh[elem[noe, 3]] };
                double del(double xi, double eta, double zeta) => z[0] * (1 - xi - eta - zeta) + z[1] * xi + z[2] * eta + z[3] * zeta
                - u(x[0] + (x[1] - x[0]) * xi + (x[2] - x[0]) * eta + (x[3] - x[0]) * zeta,
                    y[0] + (y[1] - y[0]) * xi + (y[2] - y[0]) * eta + (y[3] - y[0]) * zeta,
                    t[0] + (t[1] - t[0]) * xi + (t[2] - t[0]) * eta + (t[3] - t[0]) * zeta);
                double del2(double xi, double eta, double zeta) => del(xi, eta, zeta) * del(xi, eta, zeta);
                res += Gauss3D(del2) * jacobi;

            }
            return Math.Sqrt(res);
        }

        // Error estimate of sp fem in W(0, T) or L^2(0, T; H^1_0(Omega))
        public double ErrorW(Func<double, double, double, double > axx, Func<double, double, double, double> ayy,
            Func<double, double, double, double> dxu, Func<double, double, double, double> dyu,
            double[] uh, double[,] node, int[,] elem, double jacobi)
        {
            double res = 0;
            int NoE = elem.GetLength(0);
            for (int noe = 0; noe < NoE; noe++)
            {
                double[] x = new double[] { node[elem[noe, 0], 0], node[elem[noe, 1], 0], node[elem[noe, 2], 0], node[elem[noe, 3], 0] };
                double[] y = new double[] { node[elem[noe, 0], 1], node[elem[noe, 1], 1], node[elem[noe, 2], 1], node[elem[noe, 3], 1] };
                double[] t = new double[] { node[elem[noe, 0], 2], node[elem[noe, 1], 2], node[elem[noe, 2], 2], node[elem[noe, 3], 2] };

                double[] z = new double[] { uh[elem[noe, 0]], uh[elem[noe, 1]], uh[elem[noe, 2]], uh[elem[noe, 3]] };
                double[,] mat = new double[,] { {x[1] - x[0], x[2] - x[0], x[3] - x[0]},
                                            {y[1] - y[0], y[2] - y[0], y[3] - y[0]},
                                            {t[1] - t[0], t[2] - t[0], t[3] - t[0]}};

                //Tinh dao ham [xi, eta, zeta] theo bien x
                double D_xi_x = mat[1, 1] * mat[2, 2] - mat[2, 1] * mat[1, 2];
                double D_eta_x = -mat[1, 0] * mat[2, 2] + mat[2, 0] * mat[1, 2];
                double D_zeta_x = mat[1, 0] * mat[2, 1] - mat[2, 0] * mat[1, 1];

                //Tinh dao ham [xi, eta, zeta] theo bien y
                double D_xi_y = -mat[0, 1] * mat[2, 2] + mat[2, 1] * mat[0, 2];
                double D_eta_y = mat[0, 0] * mat[2, 2] - mat[2, 0] * mat[0, 2];
                double D_zeta_y = -mat[0, 0] * mat[2, 1] + mat[2, 0] * mat[0, 1];

                // axx * dx(uh) * dx(uh) + ayy * dy(uh) * dy(uh) => a real number
                // note: sign(J), I'm ignoring it...
                double dxuh = D_xi_x * (z[1] - z[0]) + D_eta_x * (z[2] - z[0]) + D_zeta_x * (z[3] - z[0]);
                double dyuh = D_xi_y * (z[1] - z[0]) + D_eta_y * (z[2] - z[0]) + D_zeta_y * (z[3] - z[0]);

                double x_hat(double xi, double eta, double zeta) => x[0] + (x[1] - x[0]) * xi + (x[2] - x[0]) * eta + (x[3] - x[0]) * zeta;
                double y_hat(double xi, double eta, double zeta) => y[0] + (y[1] - y[0]) * xi + (y[2] - y[0]) * eta + (y[3] - y[0]) * zeta;
                double t_hat(double xi, double eta, double zeta) => t[0] + (t[1] - t[0]) * xi + (t[2] - t[0]) * eta + (t[3] - t[0]) * zeta;

                double axx_hat(double xi, double eta, double zeta) => axx(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta), t_hat(xi, eta, zeta));
                double ayy_hat(double xi, double eta, double zeta) => ayy(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta), t_hat(xi, eta, zeta));

                double dxu_hat(double xi, double eta, double zeta) => dxu(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta), t_hat(xi, eta, zeta));
                double dyu_hat(double xi, double eta, double zeta) => dyu(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta), t_hat(xi, eta, zeta));

                double del(double xi, double eta, double zeta) =>
                    axx_hat(xi, eta, zeta) * dxu_hat(xi, eta, zeta) + ayy_hat(xi, eta, zeta) * dyu_hat(xi, eta, zeta) - (dxuh + dyuh) / jacobi;
                double del2(double xi, double eta, double zeta) => del(xi, eta, zeta) * del(xi, eta, zeta);
 
                res += Gauss3D(del2) * jacobi;

            }
            return Math.Sqrt(res);
        }
    }
}
