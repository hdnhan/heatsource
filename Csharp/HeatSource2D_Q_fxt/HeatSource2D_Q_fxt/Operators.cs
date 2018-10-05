using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HeatSource2D_Q_fxt {
    class Operators {
        // Phần định nghĩa các toán tử trong vector hay ma trận

        // Tạo vector có giá trị toàn là 1 (one)
        public double[] Ones(int LengthOfVector) {
            double[] res = new double[LengthOfVector];
            for (int i = 0; i < LengthOfVector; i++) {
                res[i] = 1;
            }
            return res;
        }

        // Tạo ma trận có giá trị toàn là 0 (zero)
        public double[,] Zeros(int row, int col) {
            double[,] res = new double[row + 1, col + 1];
            for (int i = 0; i <= row; i++) {
                for (int j = 0; j <= col; j++)
                    res[i, j] = 0;
            }
            return res;
        }

        // Phần định nghĩa các toán tử trong vector hay ma trận

        // Substract
        // z = x-y
        public double[] Sub(double[] x, double[] y) {
            double[] res = new double[x.Length];
            for (int i = 0; i < x.Length; i++) {
                res[i] = x[i] - y[i];
            }
            return res;
        }

        // C = A-B
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
        // z = x+y
        public double[] Add(double[] x, double[] y) {
            double[] res = new double[x.Length];
            for (int i = 0; i < x.Length; i++) {
                res[i] = x[i] + y[i];
            }
            return res;
        }

        // C = A+B
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


        // y = a*x
        public double[] Mul(double a, double[] x) {
            double[] res = new double[x.Length];
            for (int i = 0; i < x.Length; i++) {
                res[i] = a * x[i];
            }
            return res;
        }

        // Y = a*X
        public double[,] Mul(double a, double[,] X) {
            double[,] res = new double[X.GetLength(0), X.GetLength(1)];
            for (int i = 0; i < X.GetLength(0); i++) {
                for (int j = 0; j < X.GetLength(1); j++) {
                    res[i, j] = a * X[i, j];
                }
            }
            return res;
        }

        // y = A*x
        public double[] Mul(double[,] A, double[] x) {
            double[] res = new double[x.Length];
            for (int i = 0; i < A.GetLength(0); i++) {
                for (int j = 0; j < A.GetLength(1); j++) {
                    res[i] += A[i, j] * x[j];
                }
            }
            return res;
        }

        // a*mapA 
        public Dictionary<Tuple<int, int>, double> Mul(double a, Dictionary<Tuple<int, int>, double> mapA) {
            var res = new Dictionary<Tuple<int, int>, double>();
            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapA) {
                res.Add(KeyValue.Key, a * mapA[KeyValue.Key]);
            }
            return res;
        }

        // mapA*x
        public double[] Mul(Dictionary<Tuple<int, int>, double> mapA, double[] x) {
            double[] res = new double[x.Length];
            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapA) {
                res[KeyValue.Key.Item1] += mapA[KeyValue.Key] * x[KeyValue.Key.Item2];
            }
            return res;
        }

        // mapA'*x
        public double[] MulI(Dictionary<Tuple<int, int>, double> mapA, double[] x) {
            double[] res = new double[x.Length];
            foreach (KeyValuePair<Tuple<int, int>, double> KeyValue in mapA) {
                res[KeyValue.Key.Item2] += mapA[KeyValue.Key] * x[KeyValue.Key.Item1];
            }
            return res;
        }


        // z = x'*y
        public double InnerProduct(double[] x, double[] y) {
            double res = 0;
            for (int i = 0; i < x.Length; i++) {
                res += x[i] * y[i];
            }
            return res;
        }


        // Lấy một cột từ ma trận A
        public double[] ColumnOfMatrix(double[,] A, int index) {
            double[] res = new double[A.GetLength(0)];
            for (int i = 0; i < res.Length; i++) {
                res[i] = A[i, index];
            }
            return res;
        }

        // Lấy một hàng từ ma trận A
        public double[] RowOfMatrix(double[,] A, int index) {
            double[] res = new double[A.GetLength(1)];
            for (int i = 0; i < res.Length; i++) {
                res[i] = A[index, i];
            }
            return res;
        }



        // flip time direction
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


        // Thay đổi giá trị của một cột nào đó
        public void ChangeValueOfColumn(ref double[,] Original, int ColumnIndex, double[] InsertValue) {
            for (int i = 0; i < InsertValue.Length; i++) {
                Original[i, ColumnIndex] = InsertValue[i];
            }
        }

        // Thay đổi nghiệm từ vecto thành ma trận
        public double[,] Sol_Vec2Mat(double[] uh, int Nx, int Nt) {
            double[,] res = new double[Nx + 1, Nt + 1];
            for (int nt = 0; nt <= Nt; nt++) {
                for (int nx = 0; nx <= Nx; nx++) {
                    res[nx, nt] = uh[nx + nt * (Nx + 1)];
                }
            }
            return res;
        }

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
