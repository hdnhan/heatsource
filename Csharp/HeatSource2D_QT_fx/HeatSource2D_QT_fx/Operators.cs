using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_QT_fx {
    class Operators {
        
        public double[] Ones(int LengthOfVector) {
            double[] res = new double[LengthOfVector];
            for (int i = 0; i < LengthOfVector; i++) {
                res[i] = 1;
            }
            return res;
        }

        public double[] Zeros(int LengthOfVector) {
            double[] res = new double[LengthOfVector];
            for (int i = 0; i < LengthOfVector; i++) {
                res[i] = 0;
            }
            return res;
        }

        // Substract
        // z = x-y
        public double[] Sub(double[] x, double[] y) {
            double[] res = new double[x.Length];
            for (int i = 0; i < x.Length; i++) {
                res[i] = x[i] - y[i];
            }
            return res;
        }

        // C = A - B
        public double[,] Sub(double[,] A, double[,] B) {
            double[,] res = new double[A.GetLength(0), A.GetLength(1)];
            for (int i = 0; i < A.GetLength(0); i++) {
                for (int j = 0; j < A.GetLength(1); j++) {
                    res[i, j] = A[i, j] - B[i, j];
                }
            }
            return res;
        }

        // mapA - mapB: có thể A vs B ko cùng size, lấy A làm kết quả, B làm vòng lặp
        public Dictionary<Tuple<int, int>, double> Sub(Dictionary<Tuple<int, int>, double> mapA, Dictionary<Tuple<int, int>, double> mapB) {
            var res = new Dictionary<Tuple<int, int>, double>();

            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapA) {
                res[KeyValue.Key] = mapA[KeyValue.Key];
            }

            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapB) {
                var key = KeyValue.Key;
                double value;
                res.TryGetValue(key, out value);
                res[key] = value - mapB[key];
            }
            return res;
        }



        //Add
        // z = x + y
        public double[] Add(double[] x, double[] y) {
            double[] res = new double[x.Length];
            for (int i = 0; i < x.Length; i++) {
                res[i] = x[i] + y[i];
            }
            return res;
        }

        // C = A + B
        public double[,] Add(double[,] A, double[,] B) {
            double[,] res = new double[A.GetLength(0), A.GetLength(1)];
            for (int i = 0; i < A.GetLength(0); i++) {
                for (int j = 0; j < A.GetLength(1); j++) {
                    res[i, j] = A[i, j] + B[i, j];
                }
            }
            return res;
        }

        // mapA + mapB
        public Dictionary<Tuple<int, int>, double> Add(Dictionary<Tuple<int, int>, double> mapA, Dictionary<Tuple<int, int>, double> mapB) {
            var res = new Dictionary<Tuple<int, int>, double>();

            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapA) {
                res[KeyValue.Key] = mapA[KeyValue.Key];
            }

            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapB) {
                var key = KeyValue.Key;
                double value;
                res.TryGetValue(key, out value);
                res[key] = value + mapB[key];
            }
            return res;
        }


        // y = a * x
        public double[] Mul(double a, double[] x) {
            double[] res = new double[x.Length];
            for (int i = 0; i < x.Length; i++) {
                res[i] = a * x[i];
            }
            return res;
        }

        // Y = a * X
        public double[,] Mul(double a, double[,] X) {
            double[,] res = new double[X.GetLength(0), X.GetLength(1)];
            for (int i = 0; i < X.GetLength(0); i++) {
                for (int j = 0; j < X.GetLength(1); j++) {
                    res[i, j] = a * X[i, j];
                }
            }
            return res;
        }

        // y = A * x
        public double[] Mul(double[,] A, double[] x) {
            double[] res = new double[x.Length];
            for (int i = 0; i < A.GetLength(0); i++) {
                for (int j = 0; j < A.GetLength(1); j++) {
                    res[i] += A[i, j] * x[j];
                }
            }
            return res;
        }

        // a * mapA 
        public Dictionary<Tuple<int, int>, double> Mul(double a, Dictionary<Tuple<int, int>, double> mapA) {
            var res = new Dictionary<Tuple<int, int>, double>();
            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapA) {
                res.Add(KeyValue.Key, a * mapA[KeyValue.Key]);
            }
            return res;
        }

        // mapA * x
        public double[] Mul(Dictionary<Tuple<int, int>, double> mapA, double[] x) {
            double[] res = new double[x.Length];
            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapA) {
                res[KeyValue.Key.Item1] += mapA[KeyValue.Key] * x[KeyValue.Key.Item2];
            }
            return res;
        }

        // mapA' * x (A transpose)
        public double[] MulI(Dictionary<Tuple<int, int>, double> mapA, double[] x) {
            double[] res = new double[x.Length];
            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapA) {
                res[KeyValue.Key.Item2] += mapA[KeyValue.Key] * x[KeyValue.Key.Item1];
            }
            return res;
        }


        // z = x' * y
        public double InnerProduct(double[] x, double[] y) {
            double res = 0;
            for (int i = 0; i < x.Length; i++) {
                res += x[i] * y[i];
            }
            return res;
        }


        // flip time direction f(x, y, t)
        public double[] flip(double[] f, int Nx, int Ny, int Nt) {
            double[] res = new double[f.Length];
            for (int nt = 0; nt <= Nt / 2; nt++) {
                for (int ny = 0; ny <= Ny; ny++) {
                    for (int nx = 0; nx <= Nx; nx++) {
                        res[nt * (Nx + 1) * (Ny + 1) + ny * (Nx + 1) + nx] = f[(Nt - nt) * (Nx + 1) * (Ny + 1) + ny * (Nx + 1) + nx];
                        res[(Nt - nt) * (Nx + 1) * (Ny + 1) + ny * (Nx + 1) + nx] = f[nt * (Nx + 1) * (Ny + 1) + ny * (Nx + 1) + nx];
                    }
                }
            }
            return res;
        }

        // flip time direction f(t)
        public double[] flip(double[] f, int Nt) {
            double[] res = new double[f.Length];
            for (int nt = 0; nt <= Nt / 2; nt++) {
                res[nt] = f[Nt - nt];
                res[Nt - nt] = f[nt];
            }
            return res;
        }





        // Delete after using all of these functions
        // In ra ma trận
        public void Print(double[,] A) {
            Console.WriteLine();
            for (int i = 0; i < A.GetLength(0); i++) {
                for (int j = 0; j < A.GetLength(1); j++) {
                    Console.Write(A[i, j] + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        public void Print(int[,] A) {
            Console.WriteLine();
            for (int i = 0; i < A.GetLength(0); i++) {
                for (int j = 0; j < A.GetLength(1); j++) {
                    Console.Write(A[i, j] + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        // In ra vector
        public void Print(double[] a) {
            Console.WriteLine();
            for (int i = 0; i < a.Length; i++) {
                Console.WriteLine(a[i]);
            }
            Console.WriteLine();
        }

        // Print
        public void Print(Dictionary<Tuple<int, int>, double> map) {
            Console.WriteLine();
            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in map) {
                Console.Write(KeyValue.Value + " ");
            }
            Console.WriteLine();
        }

        public void Print(int[] a) {
            Console.WriteLine();
            for (int i = 0; i < a.Length; i++) {
                Console.WriteLine(a[i]);
            }
            Console.WriteLine();
        }
    }
}
