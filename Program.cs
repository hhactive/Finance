using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Statistics;






namespace Financial_Modeling
{
    class Program
    {
        
        public static double Pi2 = 6.28318530717959;
        public static bool store1 = false;
        public static double z = 0;
        static void Main(string[] args)
        {
            //ProjectC12();
            //ProjectC3();
            //ProjectC4();
            ProjectD();
        }
        //--------------------Project C----------------------------------------------------------------
        //---------------------------------------------------------------------------------------------
        static void ProjectC12()
        {
            Console.WriteLine("This is the Project-C 1&2.");
            int numbers = 757; //still manul input the numbers
            int StockN = 4;
            int Npaths = 20;
            double dt = 1 / 252.0;//Can be input
            double RiskFR = 0.03;  //It's a input
            double[] sigma = new double[StockN];
            double[] DivY = new double[StockN];
            double[,] Price = new double[numbers, StockN ];
            double[,] LogReturn = new double[numbers,StockN ];
            double[,] covMatrix = new double[StockN, StockN];
            double[,] Larray = new double[StockN, StockN];
            double[,] Larray2 = new double[StockN, StockN];
            double[,] Corr = new double[StockN, StockN];
            double[,,] paths = new double[numbers, Npaths, StockN];
            double[] cRV = new double[StockN];
            string FilePath = "C:\\Users\\hhact\\My Cloud\\python\\github\\C# project\\Stock.CSV";

            Price = OpenCSV(FilePath);
            //calculate logreturn and sigma
            for (int j = 0; j < StockN ; j++)
            {
                double[] tem2 = new double[numbers];
                for (int i = 0; i <numbers-1; i++)
                {
                    LogReturn[i, j] = Math.Log(Price[i, j] / Price[i + 1, j]);
                    tem2[i] = LogReturn[i, j];
                }
                sigma[j]= CalculateStd(tem2, numbers)*Math.Sqrt(252);
            }

            covMatrix =CalculateMatrixCovariance(LogReturn);
            //calculate correlation matrix
            for (int i = 0; i < StockN ; i++)
            {
                for (int j = 0; j < StockN ; j++)
                {
                    Corr[i, j] =covMatrix[i,j]*252/(sigma[i]*sigma[j]);
                    Corr[j, i] = Corr[i, j];
                }
                Corr[i, i] = 1;
            }
            Larray = Cholesky(Corr);
           
            double kk = 0;
            double p1 = 0;
            double p2 = 0;
            //Now we calculated L_matrix, we can go to next step:generate random path.
            for (int j = 0; j < Npaths; j++)
            {
                for (int i = 0; i < StockN ; i++)
                {
                    paths[0, j, i] = Price[0, i];
                }
            }

            for (int i = 1; i < numbers ; i++)
            {
                for (int j = 0; j < Npaths ; j++)
                {
                    cRV = rCGauss(Larray,StockN );
                    for (int l = 0; l < StockN ; l++)
                    {
                        p1 = (RiskFR - DivY[l] - Math.Pow(sigma[l], 2) / 2) * dt;
                        p2 = Math.Pow(dt, 0.5) * sigma[l] * cRV[l];
                        kk = Math.Exp(p1+p2);
                        paths[i, j, l] = paths[i - 1, j, l] * kk;
                            
                    }
                }
            }
            double[] AvgPrice3 = new double[4];
            for (int i = 0; i < Npaths ; i++)
            {
                AvgPrice3[0] += paths[numbers - 1,i, 0] / Npaths;
                AvgPrice3[1] += paths[numbers - 1, i, 1] / Npaths;
                AvgPrice3[2] += paths[numbers - 1, i, 2] / Npaths;
                AvgPrice3[3] += paths[numbers - 1, i, 3] / Npaths;
            }
            for (int i = 0; i < 4; i++)
            {
                Console.WriteLine("The average price of stock {0} is {1}.",i+1,AvgPrice3[i]);
            }
            //calculate option price
            Console.WriteLine("Please enter the time to maturity(Year).");
            double ttm = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Please enter the risk-free rate");
            double rf = Convert.ToDouble(Console.ReadLine());
            for (int i = 0; i < 4; i++)
            {
                double s0 = Price[0,i];

                Console.WriteLine("Please enter the strike price of stock {0}.",i+1);
                double X = Convert.ToDouble(Console.ReadLine());
                double vol = sigma[i];

                double d1 = (Math.Log(s0 / X) + (rf + (0.5 * Math.Pow(vol, 2))) * ttm) / (vol * Math.Sqrt(ttm));
                double d2 = d1 - (vol * Math.Sqrt(ttm));

                double CallOption = (s0 * CND(d1)) - (X * Math.Exp(-rf * ttm) * CND(d2));
                double Putoption = CallOption - s0 + (X * Math.Exp(-rf * ttm));
                Console.WriteLine("Call Option Price of Stock {0}:  {1}.",i+1, CallOption);
                Console.WriteLine("Put Option Price of Stock {0}:  {1}.",i+1,Putoption);
            }

        }
        
        //Part 1 for data input
        public static double[,] OpenCSV(string filePath)
        {
            System.Text.Encoding encoding = GetType(filePath);
            double[,] pp = new double[760, 4];
            System.IO.FileStream fs = new System.IO.FileStream(filePath, System.IO.FileMode.Open, System.IO.FileAccess.Read);
            System.IO.StreamReader sr = new System.IO.StreamReader(fs, encoding);

            string strLine = "";
            string[] aryLine = null;

            int Count = 0;

            while ((strLine = sr.ReadLine()) != null)
            {
                aryLine = strLine.Split(',');
                pp[Count, 0] = Convert.ToDouble(aryLine[0]);
                pp[Count, 1] = Convert.ToDouble(aryLine[1]);
                pp[Count, 2] = Convert.ToDouble(aryLine[2]);
                pp[Count, 3] = Convert.ToDouble(aryLine[3]);
                Count += 1;
            }

            sr.Close();
            fs.Close();
            return pp;
        }
        public static System.Text.Encoding GetType(string FILE_NAME)
        {
            System.IO.FileStream fs = new System.IO.FileStream(FILE_NAME, System.IO.FileMode.Open,
                System.IO.FileAccess.Read);
            System.Text.Encoding r = GetType(fs);
            fs.Close();
            return r;
        }

        public static System.Text.Encoding GetType(System.IO.FileStream fs)
        {
            byte[] Unicode = new byte[] { 0xFF, 0xFE, 0x41 };
            byte[] UnicodeBIG = new byte[] { 0xFE, 0xFF, 0x00 };
            byte[] UTF8 = new byte[] { 0xEF, 0xBB, 0xBF }; //BOM  
            System.Text.Encoding reVal = System.Text.Encoding.Default;

            System.IO.BinaryReader r = new System.IO.BinaryReader(fs, System.Text.Encoding.Default);
            int i;
            int.TryParse(fs.Length.ToString(), out i);
            byte[] ss = r.ReadBytes(i);
            if (IsUTF8Bytes(ss) || (ss[0] == 0xEF && ss[1] == 0xBB && ss[2] == 0xBF))
            {
                reVal = System.Text.Encoding.UTF8;
            }
            else if (ss[0] == 0xFE && ss[1] == 0xFF && ss[2] == 0x00)
            {
                reVal = System.Text.Encoding.BigEndianUnicode;
            }
            else if (ss[0] == 0xFF && ss[1] == 0xFE && ss[2] == 0x41)
            {
                reVal = System.Text.Encoding.Unicode;
            }
            r.Close();
            return reVal;
        }

        private static bool IsUTF8Bytes(byte[] data)
        {
            int charByteCounter = 1;    
            byte curByte;
            for (int i = 0; i < data.Length; i++)
            {
                curByte = data[i];
                if (charByteCounter == 1)
                {
                    if (curByte >= 0x80)
                    { 
                        while (((curByte <<= 1) & 0x80) != 0)
                        {
                            charByteCounter++;
                        } 
                        if (charByteCounter == 1 || charByteCounter > 6)
                        {
                            return false;
                        }
                    }
                }
                else
                { 
                    if ((curByte & 0xC0) != 0x80)
                    {
                        return false;
                    }
                    charByteCounter--;
                }
            }
            if (charByteCounter > 1)
            {
                throw new Exception("Not the correct type");
            }
            return true;
        }

        //Part 1 ending
        //part 2 matrix 
        public static double[,] CalculateMatrixCovariance(double[,] matrix)
        {

            return matrix.Covariance();



        }
        public static double[,] Cholesky(double[,] a)
        {
            int n = (int)Math.Sqrt(a.Length);

            double[,] ret = new double[n, n];
            for (int r = 0; r < n; r++)
                for (int c = 0; c <= r; c++)
                {
                    if (c == r)
                    {
                        double sum = 0;
                        for (int j = 0; j < c; j++)
                        {
                            sum += ret[c, j] * ret[c, j];
                        }
                        ret[c, c] = Math.Sqrt(a[c, c] - sum);
                    }
                    else
                    {
                        double sum = 0;
                        for (int j = 0; j < c; j++)
                            sum += ret[r, j] * ret[c, j];
                        ret[r, c] = 1.0 / ret[c, c] * (a[r, c] - sum);
                    }
                }

            return ret;
        }
        //part 2 ending
        //part 3 rCGauss
        public static double[] rCGauss(double[,] rarray,int SN)
        {
            double[] tem = new double[SN];
            double[] rNV = new double[SN];
            for (int i = 0; i < SN ; i++)
            {
                tem[i] = GetRandomNumber2();
                for (int j = 0  ; j <=i; j++)
                {
                    rNV[i] = rNV[i] + rarray[i, j] * tem[j];
                }
            }
            return rNV;
        }
        static double GetRandomNumber2()
        {
            if (store1) 
            {
                store1 = false;
                return z;
            }
            else
            {
                double r1 = Math.Log(1 - GetRandomNumber(1, 7953) / 7954.0);
                double r2 = GetRandomNumber(1, 7953) / 7954.0;
                z = Math.Sqrt(-2 *r1) * Math.Cos(Pi2 * r2);
                store1 = true;
                return z * Math.Tan(Pi2 * r2);
            }

            Console.WriteLine("/r/n");
            Console.WriteLine("/r/n");
        }
        //--------------------Project C.1 C.2 end----------------------------------------------------------

        //--------------------Project C.3----------------------------------------------------------------
        
        static void ProjectC3()
        {
            Console.WriteLine("This is the project C-3");
            Console.WriteLine("Enter the number of Rows and Columns :");
            int m = Convert.ToInt32(Console.ReadLine());
            double[][] matrixA = new double[m][];
            for (int i = 0; i < matrixA.Length; i++)
            {
                matrixA[i] = new double[matrixA.Length];
            }
            Console.WriteLine("Enter the Matrix");
            for (int i = 0; i < m; i++)
            {
                string [] aryLine = null;
                string ss="";
                ss=Console.ReadLine();
                aryLine = ss.Split(' ');
                for (int j = 0; j < m; j++)
                {
                    
                        double tem=double.Parse(aryLine[j]);
                    matrixA[i][j] = tem;
                }
            }
            


            double[][] matrixP = new double[matrixA.Length][];
            for (int i = 0; i < matrixP .Length; i++)
            {
                matrixP[i] = new double[matrixA.Length ];
            }
            double[][] matrixL = new double[matrixA.Length ][];
            for (int i = 0; i < matrixL.Length; i++)
            {
                matrixL[i] = new double[matrixA.Length];
                for (int j = 0; j < matrixA.Length; j++)
                {
                    matrixL[i][j] = 0; 
                }
            }
            double[][] matrixU = new double[matrixA.Length][];
            for (int i = 0; i < matrixU.Length; i++)
            {
                matrixU[i] = new double[matrixA.Length];
            }

            for (int i = 0; i < matrixA.Length-1 ; i++)
            {
                int position = findmax(matrixA, i);
                matrixP[i][position] = 1;
                matrixA = Mswitch(matrixA, i, position);
                matrixL = Mswitch(matrixL, i, position);
                for (int j = 1; j < matrixA.Length-i  ; j++)
                {
                    matrixL[j + i][i] = matrixA[j + i][i] / matrixA[i][i];
                    for (int k = 0; k < matrixA.Length - i; k++)
                    {
                        matrixA[j + i][k+i] = matrixA[j + i][k+i] - matrixL[j + i][i] * matrixA[i][k+i];
                            
                    }
                }
            }
            matrixP[matrixP.Length - 1][0] = 1;
            PrintMatrix(matrixA, "MatrixU");
            PrintMatrix(matrixL, "MatrixL");
            PrintMatrix(matrixP, "MatrixP");

        }
        public static double[][] Mswitch(double [][] matrix,int a, int b)
        {
         
            double temp=0;
            for (int i = 0; i < matrix.Length ; i++)
            {
                temp= matrix[b][i];
                matrix[b][i] = matrix[a][i];
                matrix[a][i] = temp;

            }
            return matrix;
        }
        public static int findmax(double[][] matrix, int num)
        {
            int p=num;
            for (int i = num+1; i < matrix.Length ; i++)
            {
                if (Math.Abs( matrix[i][num])> Math.Abs( matrix[p][num])) {
                    p = i;
                }
           
            }
            return p;
        }
        //---------------------------------------------------------------------------------------------

        //--------------------Project C.4----------------------------------------------------------------
        static void ProjectC4()
        {
            Console.WriteLine("This is project C-4.");
            int StockN = 1;
            int Numbers = 249;
            Double[] Price = new double[Numbers];
            double[][] Price1 = new double[Numbers][];
            for (int i = 0; i < Price1.Length; i++)
            {
                Price1[i] = new double[1];
            }
            double[][] X = new double[Numbers][];
            for (int i = 0; i < X.Length; i++)
            {
                X[i] = new double[2];
            }
            double[][] Beta = new double[2][];
            for (int i = 0; i < Beta.Length; i++)
            {
                Beta[i] = new double[1];
            }
            string FilePath = "C:\\Users\\hhact\\My Cloud\\python\\github\\C# project\\IBM.CSV";
            Price = OpenCSV1(FilePath, Numbers);
            for (int i = 0; i < Numbers; i++)
            {
                X[i][0] = 1;
                X[i][1] = i + 1;
                //X[i][2] = Math.Cos(i + 1);
                Price1[i][0] = Price[i];
            }
            Beta = OLS(X, Price1,Numbers );
            Console.WriteLine("For the least squares approximation: Y={0}+{1}t",Beta[0][0],Beta[1][0]);

            for (int i = 0; i < Numbers ; i++)
            {
                Price1[i][0] = Math.Log(Price1[i][0]);
            }
            Beta = OLS(X, Price1, Numbers);
            Console.WriteLine("For the exponential approximation: Y={0}*e^{1}t", Math.Exp(Beta[0][0]), Beta[1][0]);
        }

        public static double[][] OLS(double[][] x,double [][] y,int numbers)
        {
            double[][] temp;

            temp = MMult(InverseMatrix(MMult(MTrans(x), x)), MMult(MTrans(x), y));
            return temp;
        }

        private static double[][] MTrans(double[][] matrix)
        {
           
            double[][] result = new double[matrix[0].Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[matrix.Length];
            }
            //rule： b[i,j]=a[j,i]
            for (int i = 0; i < result.Length; i++)
            {
                for (int j = 0; j < result[0].Length; j++)
                {
                    result[i][j] = matrix[j][i];
                }
            }
            return result;
        }
        public static double[,] transfer(double[][] array1)
        {
            
            int a = 0;
            int b = 0;
            a = array1.Length;
            b = array1[0].Length;
            double[,] temp = new double[a,b];
            for (int i = 0; i <a; i++)
            {
                for (int j = 0; j <b; j++)
                {
                    temp[i, j] = array1[i][j];
                }
            }
            return temp;
        }

        //try
        public static double[][] InverseMatrix(double[][] matrix)
        {
            //matrix should be empty
            if (matrix == null || matrix.Length == 0)
            {
                return new double[][] { };
            }

            //matrix 
            int len = matrix.Length;
            for (int counter = 0; counter < matrix.Length; counter++)
            {
                if (matrix[counter].Length != len)
                {
                    throw new Exception("matrix should be square matrix");
                }
            }

            double dDeterminant = Determinant(matrix);
            if (Math.Abs(dDeterminant) <= 1E-6)
            {
                throw new Exception("matrix can not be inverse");
            }

            double[][] result = AdjointMatrix(matrix);

            for (int i = 0; i < matrix.Length; i++)
            {
                for (int j = 0; j < matrix.Length; j++)
                {
                    result[i][j] = result[i][j] / dDeterminant;
                }
            }

            return result;
        }


        public static double Determinant(double[][] matrix)
        {
            if (matrix.Length == 0) return 0;
            else if (matrix.Length == 1) return matrix[0][0];
            else if (matrix.Length == 2)
            {
                return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
            }


            double dSum = 0, dSign = 1;
            for (int i = 0; i < matrix.Length; i++)
            {
                double[][] matrixTemp = new double[matrix.Length - 1][];
                for (int count = 0; count < matrix.Length - 1; count++)
                {
                    matrixTemp[count] = new double[matrix.Length - 1];
                }

                for (int j = 0; j < matrixTemp.Length; j++)
                {
                    for (int k = 0; k < matrixTemp.Length; k++)
                    {
                        matrixTemp[j][k] = matrix[j + 1][k >= i ? k + 1 : k];
                    }
                }

                dSum += (matrix[0][i] * dSign * Determinant(matrixTemp));
                dSign = dSign * -1;
            }

            return dSum;
        }

        public static double[][] AdjointMatrix(double[][] matrix)
        {

            double[][] result = new double[matrix.Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[matrix[i].Length];
            }

 
            for (int i = 0; i < result.Length; i++)
            {
                for (int j = 0; j < result.Length; j++)
                {

                    double[][] temp = new double[result.Length - 1][];
                    for (int k = 0; k < result.Length - 1; k++)
                    {
                        temp[k] = new double[result[k].Length - 1];
                    }


                    for (int x = 0; x < temp.Length; x++)
                    {
                        for (int y = 0; y < temp.Length; y++)
                        {
                            temp[x][y] = matrix[x < i ? x : x + 1][y < j ? y : y + 1];
                        }
                    }


                    result[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * Determinant(temp);
                }
            }



            return result;
        }

        private static void PrintMatrix(double[][] matrix, string title = "")
        {

            if (!String.IsNullOrWhiteSpace(title))
            {
                Console.WriteLine(title);
            }


            for (int i = 0; i < matrix.Length; i++)
            {
                for (int j = 0; j < matrix[i].Length; j++)
                {
                    Console.Write(matrix[i][j].ToString("0.00") + "\t");
 
                }
                Console.WriteLine();
            }


            Console.WriteLine();
        }

        public static double[][] MInv(double[][] Array1)
        {
            int m = 0;
            int n = 0;
            double[,] Array = transfer(Array1);
            m = Array.GetLength(0);
            n = Array.GetLength(1);
            double[][] array = new double[2 * m + 1][];
            for (int k = 0; k < 2 * m + 1; k++)  
            {                  
                    array[k] = new double[2*n+1];
            }
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    array[i][j] = Array[i, j];
                }
            }

            for (int k = 0; k < m; k++)
            {
                for (int t = n; t <= 2 * n; t++)
                {
                    if ((t - k) == m)
                    {
                        array[k][t] = 1.0;
                    }
                    else
                    {
                        array[k][t] = 0;
                    }
                }
            }

            for (int k = 0; k < m; k++)
            {
                if (array[k][k] != 1)
                {
                    double bs = array[k][k];
                    array[k][k] = 1;
                    for (int p = k + 1; p < 2 * n; p++)
                    {
                        array[k][p] /= bs;
                    }
                }
                for (int q = 0; q < m; q++)
                {
                    if (q != k)
                    {
                        double bs = array[q][k];
                        for (int p = 0; p < 2 * n; p++)
                        {
                            array[q][p] -= bs * array[k][p];
                        }
                    }
                    else
                    {
                        continue;
                    }
                }
            }
            double[][] NI = new double[m][];
            for (int x = 0; x < m; x++)
            {
                NI[x] = new double[2 * n];
            }
            for (int x = 0; x < m; x++)
            {
                for (int y = n; y < 2 * n; y++)
                {
                    NI[x][y - n] = array[x][y];
                }
            }
            return NI;
        }
        //multiple matrix 
        public static double[][] MMult(double[][] matrix1, double[][] matrix2)
        {
            
            int m = matrix1.Length, n = matrix2.Length, p = matrix2[0].Length;
            double[][] result = new double[m][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[p];
            }
            //rule：c[i,j]=Sigma(k=1→n,a[i,k]*b[k,j])
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < p; j++)
                {
                    //
                    for (int k = 0; k < n; k++)
                    {
                        result[i][j] += (matrix1[i][k] * matrix2[k][j]);
                    }
                }
            }
            return result;
        }
        //--------------------Project C end----------------------------------------------------------------
        //---------------------------------------------------------------------------------------------


        //--------------------Project D----------------------------------------------------------------
        /*Introduction
        This program include mainly 4 parts.
        1.Input stock price from FilePath
        2.Using Monte Carlo Simulation to simulate stock price paths.
        3.Calculating measure of risk like VaR and volatility.
        4.Using BSM Model for put/call option pricing.
        */
        //---------------------------------------------------------------------------------------------
        static void ProjectD()
        {
            Console.WriteLine("This is the Project-D");
            int numbers = 249; //still manul input the numbers
            double[] Price = new double[numbers];

            double[] LogReturn = new double[numbers];
            double[,] MC = new double[numbers, 1000];
            double AvgPrice=0;
            string FilePath = "C:\\Users\\hhact\\My Cloud\\python\\github\\C# project\\IBM.CSV";

            Price = OpenCSV1(FilePath, numbers);
            /* check out put
            for (int i = 0; i < numbers; i++)
            {
                Console.WriteLine(Price[i]);
            }
            */


            //compute the log return 
            for (int i = 0; i < numbers - 1; i++)
            {
                LogReturn[i] = Math.Log(Price[i] / Price[i + 1]);
            }

            for (int i = 0; i < 1000; i++)
            {
                MC[0, i] = Price[0];
                for (int j = 1; j < numbers; j++)
                {
                    MC[j, i] = MC[j - 1, i] * Math.Exp(LogReturn[GetRandomNumber(1, numbers) - 1]);
                }
                AvgPrice += MC[numbers-1, i] / 1000.0;
            }
            Console.WriteLine("The prediction of stock price is {0}.",AvgPrice);
            //Until now, I got the monte carlo simulation for stock price.

            //Steps 3, measure of risk

            //3-1 Sort
            List<double> LogReturn1 = LogReturn.OrderBy(x => x).ToList();
            //3-2 Value at Risk
            Console.WriteLine("Please enter the confidence level.");
            int cl = Convert.ToInt32(Math.Ceiling((numbers - 1) * Convert.ToDouble(Console.ReadLine())));
            Console.WriteLine("The VaR value is {0}.", LogReturn1[numbers - cl]);
            //comments: Here we select Ceiling not floor to take more risk.

            //3-3 Volatility 
            double vol = CalculateStd(LogReturn, numbers) * Math.Sqrt(252);
            Console.WriteLine("The standard deviation is {0}.", vol);

            //Step 4 Option Modeling
            //4-1 Get the inputs for BSM 
            double s0 = Price[0];

            Console.WriteLine("Please enter the strike price.");
            double X = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Please enter the time to maturity(Year).");
            double ttm = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Please enter the risk-free rate");
            double rf = Convert.ToDouble(Console.ReadLine());

            double d1 = (Math.Log(s0 / X) + (rf + (0.5 * Math.Pow(vol, 2))) * ttm) / (vol * Math.Sqrt(ttm));
            double d2 = d1 - (vol * Math.Sqrt(ttm));

            double CallOption = (s0 * CND(d1)) - (X * Math.Exp(-rf * ttm) * CND(d2));
            double Putoption = CallOption - s0 + (X * Math.Exp(-rf * ttm));
            Console.WriteLine("Call Option Price:{0}.", CallOption);
            Console.WriteLine("Put Option Price:{0}.", Putoption);
        }
        //Part 1 for data input
        public static double[] OpenCSV1(string filePath, int numbers)
        {
            System.Text.Encoding encoding = GetType(filePath);
            double[] pp = new double[numbers];
            System.IO.FileStream fs = new System.IO.FileStream(filePath, System.IO.FileMode.Open, System.IO.FileAccess.Read);
            System.IO.StreamReader sr = new System.IO.StreamReader(fs, encoding);

            string strLine = "";
            string[] aryLine = null;

            int Count = 0;

            while ((strLine = sr.ReadLine()) != null)
            {
                aryLine = strLine.Split(',');
                pp[Count] = Convert.ToDouble(aryLine[1]);
                Count += 1;
            }

            sr.Close();
            fs.Close();
            return pp;
        }
  

        

        

        //Part 1 ending 

        //part 2 random number
        public static int GetRandomNumber(int min, int max)
        {
            int rtn = 0;
            Random r = new Random();
            byte[] buffer = Guid.NewGuid().ToByteArray();
            int iSeed = BitConverter.ToInt32(buffer, 0);
            r = new Random(iSeed);
            rtn = r.Next(min, max + 1);
            return rtn;
        }
        //part 2 ending

        //part 3 calculate std
        public static double CalculateStd(double[] price2,int numbers)
        {
            double stdev = 0;
            double avg = price2.Average();
            double sum = 0;
            for (int i = 0; i < numbers; i++)
            {
                sum += Math.Pow(avg - price2[i],2);
            }
            stdev = Math.Sqrt(sum / numbers);
            return stdev;
        }

        // Cumulative normal distribution ,source from NVIDIA
        static double CND(double d)
        {
            const double A1 = 0.31938153;
            const double A2 = -0.356563782;
            const double A3 = 1.781477937;
            const double A4 = -1.821255978;
            const double A5 = 1.330274429;
            const double RSQRT2PI = 0.39894228040143267793994605993438;

            double
            K = 1.0 / (1.0 + 0.2316419 * Math.Abs(d));

            double
            cnd = RSQRT2PI * Math.Exp(-0.5 * d * d) *
                  (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));

            if (d > 0)
                cnd = 1.0 - cnd;
            return cnd;
            //part 3 ending
        }
        //------------------------------Project D ending------------------------------------------
        //----------------------------------------------------------------------------------------
        
    }
    
}

