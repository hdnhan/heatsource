using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpaceTime1D
{
    class Solvers
    {
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


        public double ErrorL2(Func<double, double, double> u, double[] uh, double[,] node, int[,] elem, double jacobi)
        {
            double res = 0;
            int NoE = elem.GetLength(0);
            for (int noe = 0; noe < NoE; noe++)
            {
                double[] x = new double[] { node[elem[noe, 0], 0], node[elem[noe, 1], 0], node[elem[noe, 2], 0] };
                double[] t = new double[] { node[elem[noe, 0], 1], node[elem[noe, 1], 1], node[elem[noe, 2], 1] };
                double[] z = new double[] { uh[elem[noe, 0]], uh[elem[noe, 1]], uh[elem[noe, 2]] };
                Func<double, double, double> del = (xi, tau) => z[0] * (1 - xi - tau) + z[1] * xi + z[2] * tau
                - u(x[0] * (1 - xi - tau) + x[1] * xi + x[2] * tau, t[0] * (1 - xi - tau) + t[1] * xi + t[2] * tau);
                Func<double, double, double> del2 = (xi, tau) => del(xi, tau) * del(xi, tau);
                res += Gauss2D(del2) * jacobi;

            }
            return Math.Sqrt(res);
        }

        // Error estimate of sp fem in W(0, T) or L^2(0, T; H^1_0(Omega))
        public double ErrorW(Func<double, double, double> a, Func<double, double, double> dxu, double[] uh, double[,] node, int[,] elem, double jacobi)
        {
            double res = 0;
            int NoE = elem.GetLength(0);
            for (int noe = 0; noe < NoE; noe++)
            {
                double[] x = new double[] { node[elem[noe, 0], 0], node[elem[noe, 1], 0], node[elem[noe, 2], 0] };
                double[] t = new double[] { node[elem[noe, 0], 1], node[elem[noe, 1], 1], node[elem[noe, 2], 1] };
                double[] z = new double[] { uh[elem[noe, 0]], uh[elem[noe, 1]], uh[elem[noe, 2]] };
                // dx(uh) * dx(uh) => a real number
                double dxuh = (t[2] - t[0]) * (z[1] - z[0]) - (z[2] - z[0]) * (t[1] - t[0]); // note: sign(J), I'm ignoring it...
                Func<double, double, double> del = (xi, tau) =>
                dxu(x[0] * (1 - xi - tau) + x[1] * xi + x[2] * tau, t[0] * (1 - xi - tau) + t[1] * xi + t[2] * tau) - dxuh / jacobi;
                Func<double, double, double> a_hat = (xi, tau) =>
                a(x[0] * (1 - xi - tau) + x[1] * xi + x[2] * tau, t[0] * (1 - xi - tau) + t[1] * xi + t[2] * tau) * del(xi, tau) * del(xi, tau);
                res += Gauss2D(a_hat) * jacobi;

            }
            return Math.Sqrt(res);
        }
    }
}
