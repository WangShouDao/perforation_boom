using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace FSKprogram
{
    class Program
    {
        #region 定义数组
        static double[] tpxt = new double[1000];
        static double[] ppxt = new double[1000];
        static double[] sp100 = new double[10000];
        static double[] za = new double[1000000];
        static double[] zlpp = new double[100000];
        static double[] zluu = new double[100000];
        static double[] zeee = new double[100000];
        static double[] zree = new double[100000];
        static double[] zru2 = new double[100000];
        static double[] wwwl = new double[100000];
        static double[] whhl = new double[100000];
        static double[] wlll = new double[100000];
        static double[] zwww = new double[100000];
        static double[] zwhh = new double[100000];
        static double[] zwll = new double[100000];
        #endregion

        #region 定义变量
        static double xnn = 7.0;				//多方指数n
        static double aa = 300.0;			//常数A
        static double bb = 300.0;
        static double c0 = 1.45; 			//井液的常态声速
        static double g0 = 9.8e-6;             //重力加速度
        static double t0 = 0.0;
        static double r0 = 0.0;
        static double rou0 = 1000;			//井液的常态密度

        static int kg, ii;
        static double dr, dt, b0, hend, rend;
        static double r1, r2, r3, r4;

        static int nxxx = 0;
        static int nyyy = 0;
        static int nttt = 0;
        static double ttt;
        //static double pppp = 0;
        static int n, n1;

        static int lu0, lr0, lv0, lp0, le0, lq0, lu1, lr1, lv1, lp1, le1, lq1, ltt;

        static double t, tend, pxt;
        static double[, ,] uq;

        static int nr1, nr2, nr3, nr4;

        static double p100s;
        #endregion

        static void Main(string[] args)
        {           
            FileStream fs = new FileStream("F:\\C#\\FSK_program\\FSKprogram\\zls01.dat", FileMode.Open);
            FileStream fs1 = new FileStream("F:\\C#\\FSK_program\\FSKprogram\\fskpo6.out", FileMode.Create);
            FileStream fs2 = new FileStream("F:\\C#\\FSK_program\\FSKprogram\\zl10t.out", FileMode.Create);
            FileStream fs3 = new FileStream("F:\\C#\\FSK_program\\FSKprogram\\fskp66.out", FileMode.Create);
            FileStream fs4 = new FileStream("F:\\C#\\FSK_program\\FSKprogram\\pressure.out", FileMode.Create);
            StreamReader sr;
            StreamWriter sw;
            StreamWriter sw1;

            #region 读取数据
            sr = new StreamReader(fs);
            sw = new StreamWriter(fs1);
            

            string[] str1 = sr.ReadLine().Split(',');
            kg = int.Parse(str1[1]);
            sw.WriteLine("  kg = " + kg);              //上边界为自由面或固壁。
            dt = double.Parse(str1[3]);
            sw.WriteLine("  dt = " + dt);              //时间步长
            dr = double.Parse(str1[5]);
            sw.WriteLine("  dr = " + dr);              //空间步长 
            b0 = double.Parse(str1[7]);
            sw.WriteLine("  b0 = " + b0);              //粘性系数                
            string[] str2 = sr.ReadLine().Split(',');
            hend = double.Parse(str2[1]);
            sw.WriteLine("hend = " + hend);              //井深（现在取2700米）
            rend = double.Parse(str2[3]);
            sw.WriteLine("rend = " + rend);              //计算区域的高度
                
            string[] str3 = sr.ReadLine().Split(',');
            r1 = double.Parse(str3[1]);
            sw.WriteLine("  r1 = " + r1);
            r2 = double.Parse(str3[3]);
            sw.WriteLine("  r2 = " + r2);
            r3 = double.Parse(str3[5]);
            sw.WriteLine("  r3 = " + r3);
            r4 = double.Parse(str3[7]);
            sw.WriteLine("  r4 = " + r4);

            string str4 = sr.ReadLine();                    //空行
            string[] str5 = sr.ReadLine().Split(',');
            ii = int.Parse(str5[1]);
            string str6 = sr.ReadLine();
            n = Convert.ToInt32((rend - r0) / dr);
            n1 = n + 1;
            sw.WriteLine("n = " + n + "\t\tn1 = " + n1);
            sw.WriteLine("ii = " + ii);
            sw.WriteLine("\tt(i)\t\tpx(i)");
                
            string str7 = sr.ReadLine();
            int k = 0;
            while (str7 != null)
            {
                string[] str = str7.Split(',');
                tpxt[k] = double.Parse(str[0]);
                ppxt[k] = double.Parse(str[1]);
                sw.WriteLine("{0, 8 : F4}{1, 12 : F4}", str[0], str[1]);
                str7 = sr.ReadLine();
                k++;
            }
            #endregion

            #region 参数
            lu0 = 100;
	        lr0 = lu0 + n1;					//存放k时刻网格位置		lr0=901
	        lv0 = lr0 + n1;					//存放k时刻v值			lv0=1702
	        lp0 = lv0 + n1;				    //存放k时刻压力p值		lp0=2503
	        le0 = lp0 + n1;					//存放k时刻能量e值		le0=3304
	        lq0 = le0 + n1;					//存放k时刻粘性q值		lq0=4105
	        lu1 = lq0 + n1;					//存放k+1时刻速度u值	lu1=4906
	        lr1 = lu1 + n1;					//存放k+1时刻网格位置	lr1=5707
	        lv1 = lr1 + n1;				    //存放k+1时刻v值		lv1=6508
	        lp1 = lv1 + n1;					//存放k+1时刻压力p值	lp1=7309
	        le1 = lp1 + n1;					//存放k+1时刻能量e值	le1=8110
	        lq1 = le1 + n1;					//存放k+1时刻粘性q值	lq1=8911
            ltt = lq1 + n1;	                //存放有关的输出量		ltt=9712
            #endregion

            tend = tpxt[ii-1];
            pxt = ppxt[0];
            uq = new double[n1, 6, 2];
            if (kg == 1)
                zbegin1(n1, xnn, aa, c0, hend, pxt, r0, hend, dr, g0, rou0, za);
            else
                zbegin2(n1, xnn, aa, c0, hend, pxt, r0, hend, dr, g0, rou0, za);

            za[4] = 0;
            double roue;
            for (int i = 0; i < n; i++)
            {
                roue = za[le0 + i] / za[lv0 + i];			//密度与e的乘积
                za[4] = za[4] + roue * (za[lr0 + i + 1] - za[lr0 + i]);
            }
            sw.WriteLine("zzzz  t = 0");
            sw.WriteLine("\t\t*zu**" + "\t\t\t*zr(t)**" + "\t\t***zv**" + "\t\t\t\t***zp**" + "\t\t\t***ze***");
            sw.WriteLine();
            for (int i = 0; i < n1; i++)
            {
                if (i != (i / 1 * 1))
                    break;
                sw.Write("{0, 16 :F8}", za[lu0 + i]);
                sw.Write("{0, 16 :F8}", za[lr0 + i]);
                sw.Write("{0, 20 :F8}", 1 / za[lv0 + i]);
                sw.Write("{0, 16 :F8}", za[lp0 + i]);
                sw.WriteLine("{0, 16 :F8}", za[le0 + i]);
                sw.Flush();
            }
            int nx6 = n1 * 6;
            for (int i = 0; i < nx6; i++)
            {
                za[lu1 + i] = za[lu0 + i];
            }
            

            for (int i = 0; i < n; i++)
            {
                if ((r1 >= za[lr0 + i]) && (r1 < za[lr0 + i + 1])) nr1 = i + 1;
                if ((r2 >= za[lr0 + i]) && (r2 < za[lr0 + i + 1])) nr2 = i + 1;
                if ((r3 >= za[lr0 + i]) && (r3 < za[lr0 + i + 1])) nr3 = i + 1;
                if ((r4 >= za[lr0 + i]) && (r4 < za[lr0 + i + 1])) nr4 = i + 1;
            }
            sw.WriteLine("nr1 = \t" + nr1 + "\t\t" + za[lr0 + nr1 - 1]);
            sw.WriteLine("nr2 = \t" + nr2 + "\t\t" + za[lr0 + nr2 - 1]);
            sw.WriteLine("nr3 = \t" + nr3 + "\t\t" + za[lr0 + nr3 - 1]);
            sw.WriteLine("nr4 = \t" + nr4 + "\t\t" + za[lr0 + nr4 - 1]);
            sw.WriteLine("  n = \t" + n + "\t\t" + za[lr0 + n - 1]);
            sw.WriteLine("  n1 = \t" + n1 + "\t\t" + za[lr0 + n1 - 1]);
            sw.Flush();
            sw.Close();

            t = t0;
            p100s = sp100[0];
            ttt = 0.2 - dt / 2;
            za[0] = 0;
            za[1] = 0;
            za[2] = 0;
            za[3] = 0;
            za[5] = 0;
            za[6] = 0;
            sw1 = new StreamWriter(fs4);
            bool flag = true;
            while(flag)
            {
                t = t + dt;
                nxxx = nxxx + 1;
                if (t > tend) break;
                if(kg == 1)
                    zuxroupq1(n1,b0,xnn,aa,bb,c0,dt,dr,hend,rou0,g0,pxt,za);
                else
                    zuxroupq2(n1,b0,xnn,aa,bb,c0,dt,dr,hend,rou0,g0,pxt,za);
                za[3] = za[3] + pxt * (za[lr1] - za[lr0]);
                za[5] = rou0 * g0 * hend * 0.5 * g0 * t * t;
                za[6] = za[6] + 0.5 * (za[lp1 + n -1] + za[lp0 + n -1]) * (za[lr1 + n] - za[lr0 + n]);
                bool flag1 = true;
                while (flag1)
                {
                    if (nxxx != (nxxx / 100 * 100)) break;              
                    if (t < ttt)
                        break;
                    else if (t <= 100)
                        ttt = ttt + 0.2;
                    else
                        ttt = ttt + 1;

                    #region 记录t时刻各位置处的压力
                    t = Convert.ToInt32(t);
                    if (t / 1 == t)
                    {
                        double[] temp = new double[n1];                        
                        //sw1.WriteLine("t = {0}", Convert.ToInt32(t));
                        for (int i = 0; i < n1; i++)
                        {
                            temp[i] = za[lp1 + i];
                            sw1.Write("{0, 18:F12}", za[lp1 + i]);
                            //if ((i + 1) % 5 == 0)
                               // sw1.WriteLine();
                            sw1.Flush();
                        }
                        sw1.WriteLine();
                    }
                    #endregion


                    #region 赋值操作
                    za[ltt + nttt] = t;
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lu1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lr1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lu1 + nr1 - 1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lp1 + nr1 - 1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lu1 + nr2 - 1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lp1 + nr2 - 1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lu1 + nr3 - 1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lp1 + nr3 - 1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lu1 + nr4 - 1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = za[lp1 + nr4 - 1];
                    nttt = nttt + 1;
                    za[ltt + nttt] = pxt;
                    nttt = nttt + 1;
                    za[ltt + nttt] = p100s;
                    nttt = nttt + 1;
                    #endregion

                    #region 计算做功
                    zeee[nyyy] = za[0] - za[4] + za[5];		//内能Ei(le1)+动能Ek-内能Ei(le0)+阻止自重做功Wh
                    zree[nyyy] = za[1] - za[4];				//内能Ei(le1)-内能Ei(le0)
                    zru2[nyyy] = za[2];					//动能Ek
                    zwww[nyyy] = za[3];					//燃气做功W
                    zwhh[nyyy] = za[5];
                    nyyy = nyyy + 1;
                    #endregion

                    Console.WriteLine("t =" + t + "\t\tc0 = " + c0);
                    Console.WriteLine("zee = " + (za[0] - za[4] + za[5]) + "\t\tzwww = " + za[3]);
                    flag1 = false;
                }

                //int nx6 = n1 * 6;
                for (int i = 0; i < nx6; i++)
                {
                    za[lu0 + i] = za[lu1 + i];
                }

                pxt=pppxt(ii, p100s, pxt, t, tpxt, ppxt, sp100);
                if (t >= tend)
                    flag = false;
            }

            #region 写入fs2
            int nzzz;
            sw = new StreamWriter(fs2);
            sw.WriteLine("# 画线");
            sw.WriteLine("NAME\t线条.p");
            sw.WriteLine("RANK 1");
            sw.WriteLine("DATASET 1");
            sw.WriteLine("DIM \t" + nyyy);
            sw.WriteLine("BOUNDS -1 310 -1 1020");
            sw.WriteLine("INTERLACERD");
            sw.WriteLine("time 1.0");
            sw.WriteLine("GRID t");
            sw.Flush();

            nzzz = nyyy * 13;
            int m = 0;
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}",za[ltt + i]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR u0  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 1] * 1e3);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR zr  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 2]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR unr1  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 3] * 1e3);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR pnr1  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 4]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR unr2 ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 5] * 1e3);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR pnr2  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 6]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR unr3  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 7] * 1e3);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR pnr3  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 8]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR unr4  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 9] * 1e3);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR pnr4  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 10]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR Px(t)  ONNODE");
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 11]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR p100s  ONNODE"); 
            for (int i = 0; i < nzzz; i += 13)
            {
                sw.Write("{0, 20:F12}", za[ltt + i + 12]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }		//全为0

            sw.WriteLine("SCALAR zeee  ONNODE");
            for (int i = 0; i < nyyy; i ++)
            {
                sw.Write("{0, 20:F12}", zeee[i]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR zree  ONNODE"); 
            for (int i = 0; i < nyyy; i ++)
            {
                sw.Write("{0, 20:F12}", zree[i]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR zru2  ONNODE"); 
            for (int i = 0; i < nyyy; i ++)
            {
                sw.Write("{0, 20:F12}", zru2[i]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR zwww  ONNODE");
            for (int i = 0; i < nyyy; i ++)
            {
                sw.Write("{0, 20:F12}", zwww[i]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR zwhh  ONNODE"); 
            for (int i = 0; i < nyyy; i ++)
            {
                sw.Write("{0, 20:F12}", zwhh[i]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }

            sw.WriteLine("SCALAR zwll  ONNODE");
            for (int i = 0; i < nyyy; i++)
            {
                sw.Write("{0, 20:F12}", zwll[i]);
                m++;
                if ((m % 3) == 0)
                    sw.WriteLine();
                sw.Flush();
            }
            sw.Close();
            #endregion


            #region 写入fs3
            sw = new StreamWriter(fs3);
            sw.WriteLine();
            sw.WriteLine("\tdt = " + dt + "\t\tdr = " + dr + "\t\tn = " + xnn + "\t\tA = " + aa + "\t\tc0 = " + c0);
            sw.WriteLine("    nzz = " + nzzz / 13);
            sw.WriteLine();
            sw.WriteLine("\tt" + "\t\tR" + "\t\tP*" + "\t\tP(Mpa)" + "\t\tP(Mpa)" + "\t\tP(Mpa)" + "\t\tU(mm/ms)" + "\t\t(mm/ms)" + "\t\t(mm/ms)");
            sw.WriteLine("  (ms)" + "\t\t(m)" + "\t(Mpa)" + "\t\tz=r1" + "\t\tz=r2" + "\t\tz=r3" + "\t\tz=r1" + "\t\tz=r2" + "\t\tz=r3");
            sw.WriteLine();

            for (int i = 0; i < nzzz; i += 13)
            {
                int li = ltt + i;
                sw.Write("{0, 8 :F3}", za[li]);
                sw.Write("{0, 8 :F3}", za[li + 2]);
                sw.Write("{0, 8 :F3}", za[li + 11]);
                sw.Write("{0, 8 :F3}", za[li + 4]);
                sw.Write("{0, 8 :F3}", za[li + 6]);
                sw.Write("{0, 8 :F3}", za[li + 8]);
                sw.Write("{0, 8 :F3}", za[li + 3] * (1e3));
                sw.Write("{0, 8 :F3}", za[li + 5] * (1e3));
                sw.WriteLine("{0, 8 :F3}", za[li + 7] * (1e3));
                sw.Flush();
            }
            sw.Close();
            #endregion
        }

        public static double pppxt(int ii, double p100s, double pxt, double t, double[] tpxt, double[] ppxt, double[] sp100)
        {
            int nt = 0;
            for (int i = 0; i < ii; i++)
            {
                nt = i;
                if ((t >= tpxt[i]) && (t <= tpxt[i + 1]))
                {
                    pxt = ppxt[nt] + (t - tpxt[nt]) / (tpxt[nt + 1] - tpxt[nt]) * (ppxt[nt + 1] - ppxt[nt]);
                    p100s = sp100[nt] + (t - tpxt[nt]) / (tpxt[nt + 1] - tpxt[nt]) * (sp100[nt + 1] - sp100[nt]);
                    break;
                }
            }
            return pxt;
            
        }

        public static void zuxroupq1(int n1, double b, double xnn, double aa, double bb, double c0, double dt, double dr, double hxx, double rou0, double g, double pxt, double[] za)
        {
            int ns = 0;
            int m = 0; 
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    for (int k = 0; k < n1; k++)
                    {
                        uq[k, j, i] = za[lu0 + m];
                        m++;
                    }
                }
            }

            while (ns <= 100)
            {
                za[0] = 0;
                za[1] = 0;
                za[2] = 0;
                double pmax = 0, pq0, pqs, ph0, pz0, phs, pzs, pq10, pq20, pq1s, pq2s, roue, ux, p1, px;
                double p0 = 0;
                int imax = 0; ;

                for (int i = 0; i < n1; i++)
                {
                    if (i == 0)
                    {
                        pq0 = uq[i, 3, 0] + uq[i, 5, 0];
                        pqs = uq[i, 3, 1] + uq[i, 5, 1];
                        uq[i, 0, 1] = uq[i, 0, 0] - dt / (2 * dr * rou0) * ((pqs - pxt) + (pq0 - pxt))  - dt * g;
                        uq[i, 1, 1] = uq[i, 1, 0] + dt * (uq[i, 0, 1] + uq[i, 0, 0]) / 2;
                        continue;

                    }
                    else if (i == (n1 - 1))
                    {
                        ph0 = rou0 * g * (hxx - uq[i, 1, 0]);
                        pz0 = ph0 * Math.Pow((ph0 / aa + 1), (1 / xnn));
                        phs = rou0 * g * (hxx - uq[i, 1, 1]);
                        pzs = phs * Math.Pow((phs / aa + 1), (1 / xnn));
                        pq0 = uq[i-1, 3, 0] + uq[i-1, 5, 0];
                        pqs = uq[i-1, 3, 1] + uq[i-1, 5, 1];
                        uq[i, 0, 1] = uq[i, 0, 0] - dt / (2 * dr * rou0) * ((pzs - pqs) + (pz0 - pq0)) - dt * g;
                        uq[i, 1, 1] = uq[i, 1, 0] + dt * (uq[i, 0, 1] + uq[i, 0, 0]) / 2;
                    }
                    else
                    {
                        p0 = uq[i-1, 3, 1];
                        pq10 = uq[i-1, 3, 0] + uq[i-1, 5, 0];
                        pq20 = uq[i, 3, 0] + uq[i, 5, 0];
	                    pq1s = uq[i-1, 3, 1] + uq[i-1, 5, 1];
                        pq2s = uq[i, 3, 1] + uq[i, 5, 1];
	                    uq[i, 0, 1] = uq[i, 0, 0] - dt / (2 * dr * rou0) * ((pq2s - pq1s) + (pq20 - pq10)) - dt * g;
                        uq[i, 1, 1] = uq[i, 1, 0] + dt * (uq[i, 0, 1] + uq[i, 0, 0]) / 2;
                    }

                    uq[i-1, 2, 1] = (uq[i, 1, 1] - uq[i-1, 1, 1]) / (rou0 * dr);
	                uq[i-1, 3, 1] = c0 * c0 * (1 / uq[i-1, 2, 1] - rou0) + (xnn - 1) * uq[i-1, 4, 1] / uq[i-1, 2, 1];
	                uq[i-1, 4, 1] = uq[i-1, 4, 0] - 0.5 * (uq[i-1, 3, 1] + uq[i-1, 3, 0]) * 
                        (uq[i-1, 2, 1] - uq[i-1, 2, 0]);
	                ux = uq[i, 0, 1] - uq[i-1, 0, 1];
	                uq[i-1, 5, 1] = b * b * ux / (uq[i-1, 2, 1] + uq[i-1, 2, 0]) * (ux - Math.Abs(ux));
                    roue  = uq[i-1, 4, 1] / uq[i-1, 2, 1];
                    za[1] = za[1] + roue * (uq[i, 1, 1] - uq[i-1, 1, 1]);
                    za[2] = za[2] + (0.5 / uq[i - 1, 2, 1] * Math.Pow((0.5 * (uq[i-1, 0, 1] + uq[i, 0, 1])),
                        2)) * (uq[i, 1, 1] - uq[i-1, 1, 1]);
                    za[0] = za[1] + za[2];

                    if(i == (n1 - 1))
                        break;
                    p1 = uq[i-1, 3, 1];
                    px = Math.Abs(p1 - p0);
                    if(px > pmax)
                    {
                        pmax = px;
                        imax = i - 1;
                    }                    
                }    
                ns = ns + 1;
                if (pmax < 1e-8)
                {
                    if (ns >= 3)
                        break;
                }
            }

            m = 0;
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    for (int k = 0; k < n1; k++)
                    {
                        za[lu0 + m] = uq[k, j, i];
                        m++;
                    }
                }              
            }
        }

        public static void zbegin1(int n1, double xnn, double aa, double c0, double hx, double pxt, double z0, double zend, double dz, double g, double rou0, double[] za)
        {
            //uq = new double[n1, 6, 1];
            uq[0, 0, 0] = 0;
            uq[0, 1, 0] = z0;
            double hxx = zend - z0;
            for (int i = 0; i < n; i++)
            {
                uq[i + 1, 0, 0] = 0;
                uq[i, 2, 0] = 1 / rou0;
                uq[i, 5, 0] = 0;
            }

            int ns = 0;           
            while (ns <= 100)
            {
                double pmax = 0;
                double ph, p0, p1, px;
                int imax;
                for (int i = 0; i < n; i++)
                {
                    uq[i + 1, 1, 0] = uq[i, 1, 0] + rou0 * uq[i, 2, 0] * dz;			    //数组第二列为初始化的Z			
                    ph = rou0 * g * (hxx - (uq[i, 1, 0] + uq[i + 1, 1, 0]) / 2);		//初始压强PH
                    p0 = uq[i, 3, 0];
                    uq[i, 3, 0] = ph * Math.Pow((ph / aa + 1), (1 / xnn));     	   //数组第四列为压强P的初始条件
                    uq[i, 2, 0] = 1 / (rou0 * Math.Pow((uq[i, 3, 0] / aa + 1), (1 / xnn)));		//数组第三列为密度的倒数V	
                    uq[i, 4, 0] = (uq[i, 3, 0] - c0 * c0 * (1 / uq[i, 2, 0] - rou0)) * uq[i, 2, 0]/(xnn - 1);//数组第五列为e
                    p1 = uq[i, 3, 0];
                    px = Math.Abs(p1 - p0);
                    if (px > pmax)                           //计算最大压强，记录序号
                    {
                        pmax = px;
                        imax = i - 1;
                        //Console.WriteLine("ns = " + ns + " pmax = " + pmax + " imax = " + imax);
                    }    
                }
                //Console.WriteLine("ns = " + (ns - 1) + " pmax = " + pmax + " imax = " + imax);
                ns=ns+1;
                if (pmax <= 1e-8)
                {
                    if (ns >= 3)
                        break;
                }
                 
            }

            int k = 100;
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < n1; j++)
                {
                    za[k] = uq[j, i, 0];                    
                    k++;
                }
            }
            
        }

        public static void zbegin2(int n1, double xnn, double aa, double c0, double hx, double pxt, double z0, double zend, double dz, double g, double rou0, double[] za)
        {
            //uq = new double[n1, 6, 1];
            uq[0, 0, 0] = 0;
            uq[0, 1, 0] = z0;
            double hxx = zend - z0;
            for (int i = 0; i < n; i++)
            {
                uq[i + 1, 0, 0] = 0;
                uq[i, 2, 0] = 1 / rou0;
                uq[i, 5, 0] = 0;
            }

            int ns = 0;
            bool Flag = true;
            while (Flag)
            {
                double pmax = 0;
                double ph, p0, p1, px;
                int imax = 0;
                for (int i = 0; i < n; i++)
                {
                    uq[i + 1, 1, 0] = uq[i, 1, 0] + rou0 * uq[i, 2, 0] * dz;			    //数组第二列为初始化的Z			
                    ph = rou0 * g * (hxx - (uq[i, 1, 0] + uq[i + 1, 1, 0]) / 2);		//初始压强PH
                    p0 = uq[i, 3, 0];
                    uq[i, 3, 0] = ph * Math.Pow((ph / aa + 1), (1 / xnn));     	   //数组第四列为压强P的初始条件
                    uq[i, 2, 0] = 1 / (rou0 * Math.Pow((uq[i, 3, 0] / aa + 1), (1 / xnn)));		//数组第三列为密度的倒数V	
                    uq[i, 4, 0] = (uq[i, 3, 0] - c0 * c0 * (1 / uq[i, 2, 0] - rou0)) * uq[i, 2, 0] / (xnn - 1);//数组第五列为e
                    p1 = uq[i, 3, 0];
                    px = Math.Abs(p1 - p0);
                    if (px > pmax)                           //计算最大压强，记录序号
                    {
                        pmax = px;
                        imax = i - 1;
                        Console.WriteLine("ns = " + ns + " pmax = " + pmax + " imax = " + imax);
                    }
                }
                
                ns = ns + 1;
                if (ns < 100)
                {
                    if (pmax <= 1e-8)
                        if (ns >= 3)
                            Flag = false;
                }
                else
                {
                    Flag = false;
                    Console.WriteLine("ns = " + (ns - 1) + " pmax = " + pmax + " imax = " + imax);
                }

            } 

            int k = 100;
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < n1; j++)
                {
                    za[k++] = uq[j, i, 0];
                }
            }
        }

        public static void zuxroupq2(int n1, double b, double xnn, double aa, double bb, double c0, double dt, double dr, double hxx, double rou0, double g, double pxt, double[] za)
        {
        }



    }
}
