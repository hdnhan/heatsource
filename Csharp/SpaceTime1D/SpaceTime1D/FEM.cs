using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpaceTime1D
{
    class FEM
    {
        public double[] SolveFEM(double[,] node, int[,] elem, int[] dirichlet, double jacobi,
            Func<double, double, double> a, Func<double, double> dxu0, Func<double, double, double> F)
        {
            Operators oprs = new Operators();
            Solvers slvs = new Solvers();

            int NoN = node.GetLength(0); // Number of Nodes
            int NoE = elem.GetLength(0); // Number of Elements
            double[] uh = new double[NoN]; // approximate solution

            // A * uh = b  
            Dictionary<Tuple<int, int>, double> A = new Dictionary<Tuple<int, int>, double>();
            double[] b = new double[NoN]; // Left hand side of the equation
            LinearEquation(A, ref b, node, elem, jacobi, a, dxu0, F);

            // Remove nodes on dirichlet boundary
            for (int i = 0; i < dirichlet.Length; i++)
            {
                b[dirichlet[i]] = 0;
            }

            // Remove nodes on dirichlet boundary
            for (int i = 0; i < dirichlet.Length; i++)
            {
                A.ToList().Where(pair => pair.Key.Item1 == dirichlet[i]).ToList().ForEach(pair => A.Remove(pair.Key));
                A.ToList().Where(pair => pair.Key.Item2 == dirichlet[i]).ToList().ForEach(pair => A.Remove(pair.Key));
            }
            uh = slvs.BiCGM(A, b);
            return uh;
        }

        public void LinearEquation(Dictionary<Tuple<int, int>, double> A, ref double[] b, double[,] node, int[,] elem, double jacobi,
            Func<double, double, double> a, Func<double, double> dxu0, Func<double, double, double> F)
        {
            //Operators ops = new Operators();
            int NoE = elem.GetLength(0);
            for (int noe = 0; noe < NoE; noe++)
            {
                double[] x = new double[] { node[elem[noe, 0], 0], node[elem[noe, 1], 0], node[elem[noe, 2], 0] };
                double[] t = new double[] { node[elem[noe, 0], 1], node[elem[noe, 1], 1], node[elem[noe, 2], 1] };
                double[,] loc = local(x, t, jacobi, a, dxu0, F);
                for (int i = 0; i < 3; i++)
                {
                    b[elem[noe, i]] += loc[3, i];
                    for (int j = 0; j < 3; j++)
                    {
                        var key = new Tuple<int, int>(elem[noe, i], elem[noe, j]);
                        double value;
                        A.TryGetValue(key, out value);
                        A[key] = value + loc[i, j];
                    }
                }
            }
        }


        public double[,] local(double[] x, double[] t, double jacobi,
            Func<double, double, double> a, Func<double, double> dxu0, Func<double, double, double> F)
        {
            double[,] res = new double[4, 3];
            Solvers slvs = new Solvers();

            double[,] grad_phi = new double[,] { {t[1] - t[2], x[2] - x[1]},
                                                 {t[2] - t[0], x[0] - x[2]},
                                                 {t[0] - t[1], x[1] - x[0]}};

            double sgn = Math.Sign(grad_phi[1, 0] * grad_phi[2, 1] - grad_phi[1, 1] * grad_phi[2, 0]);
            //Console.WriteLine(sgn);
            Func<double, double, double> a_hat = (xi, tau) 
                => a(x[0] * (1 - xi - tau) + x[1] * xi + x[2] * tau, 
                     t[0] * (1 - xi - tau) + t[1] * xi + t[2] * tau);

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    res[i, j] = sgn * grad_phi[j, 1] / 6 + slvs.Gauss2D(a_hat) * grad_phi[i, 0] * grad_phi[j, 0] / (jacobi);
                }
            }

            Func<double, double, double> dxu0_hat = (xi, tau) 
                => dxu0(x[0] * (1 - xi - tau) + x[1] * xi + x[2] * tau);
            Func<double, double, double> F_hat = (xi, tau)
                => F(x[0] * (1 - xi - tau) + x[1] * xi + x[2] * tau,
                     t[0] * (1 - xi - tau) + t[1] * xi + t[2] * tau);

            Func<double, double, double> F_hat1 = (xi, tau) 
                => F_hat(xi, tau) * (1 - xi - tau) - sgn * a_hat(xi, tau) * dxu0_hat(xi, tau) * grad_phi[0, 0] / jacobi;
            Func<double, double, double> F_hat2 = (xi, tau) 
                => F_hat(xi, tau) * xi - sgn * a_hat(xi, tau) * dxu0_hat(xi, tau) * grad_phi[1, 0] / jacobi;
            Func<double, double, double> F_hat3 = (xi, tau) 
                => F_hat(xi, tau) * tau - sgn * a_hat(xi, tau) * dxu0_hat(xi, tau) * grad_phi[2, 0] / jacobi;

            res[3, 0] = slvs.Gauss2D(F_hat1) * jacobi;
            res[3, 1] = slvs.Gauss2D(F_hat2) * jacobi;
            res[3, 2] = slvs.Gauss2D(F_hat3) * jacobi;
            return res;
        }



        public void Domain(ref double[,] node, ref int[,] elem, ref int[] dirichlet, double[] xlim, double T, int Nx, int Nt)
        {
            double hx = (xlim[1] - xlim[0]) / Nx;
            double ht = T / Nt; // time step

            // nodes
            int i = 0;
            int diri = 0;
            for (int nt = 0; nt <= Nt; nt++)
            {
                for (int nx = 0; nx <= Nx; nx++)
                {
                    node[i, 0] = xlim[0] + nx * hx;
                    node[i, 1] = nt * ht;
                    if (nx == 0 || nx == Nx || nt == 0)
                    {
                        node[i, 2] = 1;
                        dirichlet[diri] = i;
                        diri++;
                    }
                    else
                        node[i, 2] = 0;
                    i++;
                }
            }
            //Console.WriteLine(i);
            // elements
            i = 0;
            for (int nt = 0; nt < Nt; nt++)
            {
                for (int nx = 0; nx < Nx; nx++)
                {
                    /*if ((nx == 0 && nt == Nt - 1) || (nt == 0 && nx == Nx - 1))
                    {
                        elem[i, 0] = nx + nt * (Nx + 1) + 1;
                        elem[i, 1] = nx + nt * (Nx + 1) + Nx + 2;
                        elem[i, 2] = nx + nt * (Nx + 1) + Nx + 1;


                        elem[i + 1, 0] = nx + nt * (Nx + 1) + Nx + 1;
                        elem[i + 1, 1] = nx + nt * (Nx + 1);
                        elem[i + 1, 2] = nx + nt * (Nx + 1) + 1;
                    }
                    else
                    {*/
                    elem[i, 0] = nx + nt * (Nx + 1);
                    elem[i, 1] = nx + nt * (Nx + 1) + 1;
                    elem[i, 2] = nx + nt * (Nx + 1) + Nx + 2;


                    elem[i + 1, 0] = nx + nt * (Nx + 1) + Nx + 2;
                    elem[i + 1, 1] = nx + nt * (Nx + 1) + Nx + 1;
                    elem[i + 1, 2] = nx + nt * (Nx + 1);
                    //}
                    //Console.WriteLine(i + ": " + elem[i, 0] + "," + elem[i, 1] + "," + elem[i, 2]);
                    //Console.WriteLine(i + 1 + ": " + elem[i + 1, 0] + "," + elem[i + 1, 1] + "," + elem[i + 1, 2]);
                    i += 2;
                }
            }
        }
    }
}
