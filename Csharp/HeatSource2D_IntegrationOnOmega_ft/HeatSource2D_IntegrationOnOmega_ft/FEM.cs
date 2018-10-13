using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_IntegrationOnOmega_ft {
    class FEM {
        public double[] SolveFEM(double[,] node, int[,] elem, int[] dirichlet, int Nx, int Ny, double jacobi,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy,
            Func<double, double, double> dxu0, Func<double, double, double> dyu0,
            double[] f, Func<double, double, double, double> q) {

            Operators oprs = new Operators();
            Solvers slvs = new Solvers();

            int NoN = node.GetLength(0); // Number of Nodes
            int NoE = elem.GetLength(0); // Number of Elements
            double[] uh = new double[NoN]; // approximate solution

            // A * uh = b  
            Dictionary<Tuple<int, int>, double> A = new Dictionary<Tuple<int, int>, double>();
            double[] b = new double[NoN]; // Left hand side of the equation
            LinearEquation(A, ref b, node, elem, Nx, Ny, jacobi, axx, ayy, dxu0, dyu0, f, q);

            // Remove nodes on dirichlet boundary
            for (int i = 0; i < dirichlet.Length; i++) {
                b[dirichlet[i]] = 0;
            }

            // Remove nodes on dirichlet boundary
            for (int i = 0; i < dirichlet.Length; i++) {
                A.ToList().Where(pair => pair.Key.Item1 == dirichlet[i]).ToList().ForEach(pair => A.Remove(pair.Key));
                A.ToList().Where(pair => pair.Key.Item2 == dirichlet[i]).ToList().ForEach(pair => A.Remove(pair.Key));
            }
            uh = slvs.BiCGM(A, b);
            return uh;
        }


        public void LinearEquation(Dictionary<Tuple<int, int>, double> A, ref double[] b, double[,] node, int[,] elem, int Nx, int Ny, double jacobi,
           Func<double, double, double, double> axx, Func<double, double, double, double> ayy,
           Func<double, double, double> dxu0, Func<double, double, double> dyu0,
           double[] f, Func<double, double, double, double> q) {

            int NoE = elem.GetLength(0);
            for (int noe = 0; noe < NoE; noe++) {
                double[] x = new double[] { node[elem[noe, 0], 0], node[elem[noe, 1], 0], node[elem[noe, 2], 0], node[elem[noe, 3], 0] };
                double[] y = new double[] { node[elem[noe, 0], 1], node[elem[noe, 1], 1], node[elem[noe, 2], 1], node[elem[noe, 3], 1] };
                double[] t = new double[] { node[elem[noe, 0], 2], node[elem[noe, 1], 2], node[elem[noe, 2], 2], node[elem[noe, 3], 2] };
                double[] ff = new double[] { f[elem[noe, 0] / ((Nx + 1) * (Ny + 1))], f[elem[noe, 1] / ((Nx + 1) * (Ny + 1))],
                    f[elem[noe, 2] / ((Nx + 1) * (Ny + 1))], f[elem[noe, 3] / ((Nx + 1) * (Ny + 1))]};
                double[,] loc = local(x, y, t, jacobi, axx, ayy, dxu0, dyu0, ff, q);
                for (int i = 0; i < 4; i++) {
                    b[elem[noe, i]] += loc[4, i];
                    for (int j = 0; j < 4; j++) {
                        var key = new Tuple<int, int>(elem[noe, i], elem[noe, j]);
                        double value;
                        A.TryGetValue(key, out value);
                        A[key] = value + loc[i, j];
                    }
                }
            }
        }


        public double[,] local(double[] x, double[] y, double[] t, double jacobi,
            Func<double, double, double, double> axx, Func<double, double, double, double> ayy,
            Func<double, double, double> dxu0, Func<double, double, double> dyu0,
            double[] f, Func<double, double, double, double> q) {

            double[,] res = new double[5, 4];
            Solvers slvs = new Solvers();

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

            //Tinh dao ham [xi, eta, zeta] theo bien t
            double D_xi_z = mat[0, 1] * mat[1, 2] - mat[1, 1] * mat[0, 2];
            double D_eta_z = -mat[0, 0] * mat[1, 2] + mat[1, 0] * mat[0, 2];
            double D_zeta_z = mat[0, 0] * mat[1, 1] - mat[1, 0] * mat[0, 1];

            double[,] grad_varphi = new double[,]{{-D_xi_x-D_eta_x-D_zeta_x, D_xi_x, D_eta_x, D_zeta_x},
                                                  {-D_xi_y-D_eta_y-D_zeta_y, D_xi_y, D_eta_y, D_zeta_y},
                                                  {-D_xi_z-D_eta_z-D_zeta_z, D_xi_z, D_eta_z, D_zeta_z}};

            double Jacobi = mat[0, 0] * D_xi_x + mat[0, 1] * D_eta_x + mat[0, 2] * D_zeta_x;
            double sgn = Math.Sign(Jacobi);
            //Console.WriteLine(sgn);

            double x_hat(double xi, double eta, double zeta) => x[0] + (x[1] - x[0]) * xi + (x[2] - x[0]) * eta + (x[3] - x[0]) * zeta;
            double y_hat(double xi, double eta, double zeta) => y[0] + (y[1] - y[0]) * xi + (y[2] - y[0]) * eta + (y[3] - y[0]) * zeta;
            double t_hat(double xi, double eta, double zeta) => t[0] + (t[1] - t[0]) * xi + (t[2] - t[0]) * eta + (t[3] - t[0]) * zeta;

            double axx_hat(double xi, double eta, double zeta) => axx(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta), t_hat(xi, eta, zeta));
            double ayy_hat(double xi, double eta, double zeta) => ayy(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta), t_hat(xi, eta, zeta));

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    res[i, j] = sgn * grad_varphi[2, j] / 24 + (slvs.Gauss3D(axx_hat) * grad_varphi[0, i] * grad_varphi[0, j] +
                                                                slvs.Gauss3D(ayy_hat) * grad_varphi[1, i] * grad_varphi[1, j]) / jacobi;
                }
            }

            double fq_hat(double xi, double eta, double zeta) => f[0] * q(x[0], y[0], t[0]) * (1 - xi - eta - zeta) +
                f[1] * q(x[1], y[1], t[1]) * xi + f[2] * q(x[2], y[2], t[2]) * eta + f[3] * q(x[3], y[3], t[3]) * zeta;
            //double fq_hat(double xi, double eta, double zeta) => q(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta), t_hat(xi, eta, zeta));
            double dxu0_hat(double xi, double eta, double zeta) => axx_hat(xi, eta, zeta) * dxu0(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta));
            double dyu0_hat(double xi, double eta, double zeta) => ayy_hat(xi, eta, zeta) * dyu0(x_hat(xi, eta, zeta), y_hat(xi, eta, zeta));

            double F0(double xi, double eta, double zeta) => fq_hat(xi, eta, zeta) * (1 - xi - eta - zeta)
                - sgn * dxu0_hat(xi, eta, zeta) * grad_varphi[0, 0] / jacobi
                - sgn * dyu0_hat(xi, eta, zeta) * grad_varphi[1, 0] / jacobi;
            double F1(double xi, double eta, double zeta) => fq_hat(xi, eta, zeta) * xi
                - sgn * dxu0_hat(xi, eta, zeta) * grad_varphi[0, 1] / jacobi
                - sgn * dyu0_hat(xi, eta, zeta) * grad_varphi[1, 1] / jacobi;
            double F2(double xi, double eta, double zeta) => fq_hat(xi, eta, zeta) * eta
                - sgn * dxu0_hat(xi, eta, zeta) * grad_varphi[0, 2] / jacobi
                - sgn * dyu0_hat(xi, eta, zeta) * grad_varphi[1, 2] / jacobi;
            double F3(double xi, double eta, double zeta) => fq_hat(xi, eta, zeta) * zeta
                - sgn * dxu0_hat(xi, eta, zeta) * grad_varphi[0, 3] / jacobi
                - sgn * dyu0_hat(xi, eta, zeta) * grad_varphi[1, 3] / jacobi;

            res[4, 0] = slvs.Gauss3D(F0) * jacobi;
            res[4, 1] = slvs.Gauss3D(F1) * jacobi;
            res[4, 2] = slvs.Gauss3D(F2) * jacobi;
            res[4, 3] = slvs.Gauss3D(F3) * jacobi;

            return res;
        }



        public void Domain(ref double[,] node, ref int[,] elem, ref int[] dirichlet,
            double[] xlim, double[] ylim, double T, int Nx, int Ny, int Nt) {
            double hx = (xlim[1] - xlim[0]) / Nx;
            double hy = (ylim[1] - ylim[0]) / Ny;
            double ht = T / Nt; // time step

            // nodes
            int i = 0;
            int diri = 0;
            for (int nt = 0; nt <= Nt; nt++) {
                for (int ny = 0; ny <= Ny; ny++) {
                    for (int nx = 0; nx <= Nx; nx++) {
                        node[i, 0] = xlim[0] + nx * hx;
                        node[i, 1] = ylim[0] + ny * hy;
                        node[i, 2] = nt * ht;
                        if (nx == 0 || nx == Nx || ny == 0 || ny == Ny || nt == 0) {
                            dirichlet[diri] = i;
                            diri++;
                        }
                        i++;
                    }
                }
            }
            //Console.WriteLine(i);
            // elements
            i = 0;
            int[] col1 = new int[2];
            int[] col2 = new int[2];
            int[] col3 = new int[2];
            int[] col4 = new int[2];
            for (int nt = 0; nt < Nt; nt++) {
                for (int ny = 0; ny < Ny; ny++) {
                    for (int nx = 0; nx < Nx; nx++) {
                        col1[0] = nx + ny * (Nx + 1) + nt * (Nx + 1) * (Ny + 1);
                        col1[1] = col1[0] + (Nx + 1) * (Ny + 1);

                        col2[0] = col1[0] + 1;
                        col2[1] = col1[1] + 1;

                        col3[0] = col1[0] + Nx + 1;
                        col3[1] = col1[1] + Nx + 1;

                        col4[0] = col3[0] + 1;
                        col4[1] = col3[1] + 1;

                        // elem 1
                        elem[i, 0] = col1[1];
                        elem[i, 1] = col1[0];
                        elem[i, 2] = col3[1];
                        elem[i, 3] = col4[1];
                        i++;
                        // elem 2
                        elem[i, 0] = col1[0];
                        elem[i, 1] = col1[1];
                        elem[i, 2] = col2[1];
                        elem[i, 3] = col4[1];
                        i++;
                        // elem 3
                        elem[i, 0] = col2[0];
                        elem[i, 1] = col1[0];
                        elem[i, 2] = col2[1];
                        elem[i, 3] = col4[1];
                        i++;
                        // elem 4
                        elem[i, 0] = col1[0];
                        elem[i, 1] = col2[0];
                        elem[i, 2] = col4[0];
                        elem[i, 3] = col4[1];
                        i++;

                        // elem 5
                        elem[i, 0] = col3[0];
                        elem[i, 1] = col1[0];
                        elem[i, 2] = col4[0];
                        elem[i, 3] = col4[1];
                        i++;
                        // elem 6
                        elem[i, 0] = col1[0];
                        elem[i, 1] = col3[0];
                        elem[i, 2] = col3[1];
                        elem[i, 3] = col4[1];
                        i++;

                        //Console.WriteLine(i + ": " + elem[i, 0] + "," + elem[i, 1] + "," + elem[i, 2]);
                        //Console.WriteLine(i + 1 + ": " + elem[i + 1, 0] + "," + elem[i + 1, 1] + "," + elem[i + 1, 2]);
                    }
                }
            }
        }
    }
}
