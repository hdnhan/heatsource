using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpaceTime1D
{
    class Program
    {
        static void Main(string[] args)
        {
            // DATA
            Func<double, double, double> a = (x, t) => 1;

            // Ex 1:
            Func<double, double, double> u = (x, t) => Math.Cos(Math.PI * t) * Math.Sin(Math.PI * x);
            Func<double, double, double> F = (x, t) => -Math.PI * Math.Sin(Math.PI * t) * Math.Sin(Math.PI * x) + Math.PI * Math.PI * u(x, t);
            Func<double, double, double> dxu = (x, t) => Math.PI * Math.Cos(Math.PI * t) * Math.Cos(Math.PI * x);

            // Ex 2:
            //Func<double, double, double> u = (x, t) => Math.Exp(t) * Math.Sin(Math.PI * x);
            //Func<double, double, double> F = (x, t) => (1 + Math.PI * Math.PI) * u(x, t);
            //Func<double, double, double> dxu = (x, t) => Math.PI * Math.Exp(t) * Math.Cos(Math.PI * x);


            //Func<double, double, double> uD = (x, t) => 0; //have to, unuseful
            Func<double, double> u0 = (x) => u(x, 0);
            Func<double, double> dxu0 = (x) => Math.PI * Math.Cos(Math.PI * x);

            int nn = 4;
            int Nx = nn;
            int Nt = nn;
            double[] xlim = new double[] { 0, 1 };
            double T = 1;
            double hx = (xlim[1] - xlim[0]) / Nx;
            double ht = T / Nt; // time step
            double jacobi = hx * ht; // Jacobi = 2 * area(elem)

            // Domain 2D
            double[,] node = new double[(Nx + 1) * (Nt + 1), 3];
            int[,] elem = new int[2 * Nx * Nt, 3];
            int[] dirichlet = new int[Nx + 2 * Nt + 1]; // Dirichlet boundary


            Operators oprs = new Operators();
            Solvers slvs = new Solvers();
            FEM sfem = new FEM();

            sfem.Domain(ref node, ref elem, ref dirichlet, xlim, T, Nx, Nt);
            double[] uh = sfem.SolveFEM(node, elem, dirichlet, jacobi, a, dxu0, F);
            double[] del = new double[node.GetLength(0)];
            for (int i = 0; i < node.GetLength(0); i++)
            {
                uh[i] = uh[i] + u0(node[i, 0]);
                del[i] = u(node[i, 0], node[i, 1]) - uh[i];
            }

            Console.WriteLine("Error in L2: " + slvs.ErrorL2(u, uh, node, elem, jacobi));
            Console.WriteLine("Error in W: " + slvs.ErrorW(a, dxu, uh, node, elem, jacobi));
            Console.WriteLine("min and max: " + del.Min() + ", " + del.Max());
            Console.WriteLine("Done!");
            Console.ReadLine();
        }
    }
}
