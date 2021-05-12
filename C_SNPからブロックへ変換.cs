using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SELDLA
{
    class C_SNPからブロックへ変換
    {
        public C_SNPからブロックへ変換(){
            initialize();
        }

        public string V_一行前のコンティグ名;
        public int V_コンティグごとのSNP数;
        public List<int[]> I_LIST_コンティグごとのSNPの個人のジェノタイプ一覧; //基本的に-1,0,1のみのはず
        public List<int> I_LIST_コンティグ中のPOSの一覧;
        public int V_optionで指定されるブロックの長さ_bp = 10000;
        public Encoding Enc = Encoding.GetEncoding("UTF-8");
        public StreamWriter writer;
        //public string inputvcf;
        //public string output;
        //public double opt_balance=0.1;
        public double V_optionで指定されるブロック内のSNPの一致度 = 0.9; //空の値は除外して計算
        public struct distance
        {
            public int order;
            public double match_rate;
        }

        /// <summary>
        /// 空の値を除いて一致する割合を計算
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="rateOfNotNA"></param>
        /// <returns></returns>
        public distance get_dist(int[] x, int[] y, double rateOfNotNA)
        {
            int n = 0;
            int nm = 0;
            int nn = 0;
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] != -1 && y[i] != -1)
                {
                    n++;
                    if (x[i] == y[i]) { nm++; } else { nn++; };
                }
            }
            distance result = new distance();
            if (n == 0 || n/ (double) x.Length < rateOfNotNA)
            {
                result.order = 0;
                result.match_rate = 0;
            }
            else
            {
                if (nm >= nn)
                {
                    result.order = 1;
                    result.match_rate = nm / (double)n;
                }
                else
                {
                    result.order = -1;
                    result.match_rate = nn / (double)n;
                }
            }
            return result;
        }

        public void initialize()
        {
            V_一行前のコンティグ名 = "";
            V_コンティグごとのSNP数 = 0;
            I_LIST_コンティグごとのSNPの個人のジェノタイプ一覧 = new List<int[]>();
            I_LIST_コンティグ中のPOSの一覧 = new List<int>();
        }

        /// <summary>
        /// ジェノタイプファイルを読み取りブロックでまとめた情報をファイルで出力する
        /// </summary>
        /// <param name="inputtxt"></param>
        /// <param name="opt_nc"></param>
        /// <param name="opt_r"></param>
        /// <param name="rateOfNotNA"></param>
        public void run(string inputtxt, double opt_nc, int opt_r, double rateOfNotNA)
        {
            V_optionで指定されるブロック内のSNPの一致度=opt_nc;
            V_optionで指定されるブロックの長さ_bp=opt_r;
            string output=inputtxt+".ld";
            int 行数 = 0;
            string line;

            // Read the file and display it line by line.  
            System.IO.StreamReader file = new System.IO.StreamReader(inputtxt);
            while ((line = file.ReadLine()) != null)
            {
                行数++;
                //System.Console.WriteLine(line);
                if (行数 == 1)
                {
                    string[] header = line.Split("\t");
                    //Console.WriteLine(string.Join(", ", header));
                    writer = new StreamWriter(output);
                    writer.Write("#info");
                    for (int i = 0; i < header.Length; i++)
                    {
                        writer.Write("\t" + header[i]);
                    }
                    writer.WriteLine("");
                }
                else
                {
                    string[] values = line.Split("\t");
                    int num_people = values.Length - 2; //SNPファイルの#chr  posの2列分を引く
                    if (V_一行前のコンティグ名 != values[0] && V_一行前のコンティグ名 != "" && V_コンティグごとのSNP数 > 0)
                    {
                        Console.WriteLine("making blocks of " + V_一行前のコンティグ名);
                        F_コンティグごとのブロック構築(V_一行前のコンティグ名, rateOfNotNA);
                        initialize();
                    }
                    V_一行前のコンティグ名 = values[0];
                    V_コンティグごとのSNP数++;
                    int[] tempdata = new int[num_people];
                    for (int i = 2; i < values.Length; i++)
                    {
                        tempdata[i - 2] = Int32.Parse(values[i]);
                    }
                    //コンティグごとの全データをため込むリスト
                    I_LIST_コンティグごとのSNPの個人のジェノタイプ一覧.Add(tempdata);
                    I_LIST_コンティグ中のPOSの一覧.Add(Int32.Parse(values[1]));
                }
            }
            if (V_コンティグごとのSNP数 > 0)
            {
                Console.WriteLine("making blocks of " + V_一行前のコンティグ名);
                F_コンティグごとのブロック構築(V_一行前のコンティグ名, rateOfNotNA);
            }

            file.Close();
            writer.Close();
        }

        public void F_コンティグごとのブロック構築(string contig, double rateOfNotNA)
        {
            int start_pos = 0;
            List<int[]> tempdata = new List<int[]>();
            List<int> temppos = new List<int>();
            for (int i = 0; i < V_コンティグごとのSNP数; i++)
            {
                if (start_pos == 0) { start_pos = I_LIST_コンティグ中のPOSの一覧[i]; }
                if (I_LIST_コンティグ中のPOSの一覧[i] <= start_pos + V_optionで指定されるブロックの長さ_bp)
                {
                }
                else
                {
                    範囲内にあるSNPをブロックとして纏める(tempdata, temppos, contig, rateOfNotNA);
                    tempdata = new List<int[]>();
                    temppos = new List<int>();
                    start_pos = I_LIST_コンティグ中のPOSの一覧[i];
                }
                tempdata.Add(I_LIST_コンティグごとのSNPの個人のジェノタイプ一覧[i]);
                temppos.Add(I_LIST_コンティグ中のPOSの一覧[i]);
            }
            範囲内にあるSNPをブロックとして纏める(tempdata, temppos, contig, rateOfNotNA);
        }

        public void 範囲内にあるSNPをブロックとして纏める(List<int[]> I_SNPごとのジェノタイプ一覧, List<int> I_POS一覧, string contig, double rateOfNotNA)
        {
            int V_SNP数 = I_POS一覧.Count;
            double[,] distmatrix = new double[V_SNP数, V_SNP数];
            int[,] corder = new int[V_SNP数, V_SNP数];
            bool[,] flagclust = new bool[V_SNP数, V_SNP数];
            bool[] flagid = new bool[V_SNP数];
            //全SNP間の距離を計算
            for (int i = 0; i < V_SNP数; i++)
            {
                for (int j = 0; j < V_SNP数; j++)
                {
                    // i==jの場合は0で初期化されている
                    if (i > j)
                    {
                        distance tempdist = get_dist(I_SNPごとのジェノタイプ一覧[i], I_SNPごとのジェノタイプ一覧[j], rateOfNotNA);
                        distmatrix[i, j] = tempdist.match_rate;
                        distmatrix[j, i] = tempdist.match_rate;
                        corder[i, j] = tempdist.order;
                        corder[j, i] = tempdist.order;
                        if (tempdist.match_rate >= V_optionで指定されるブロック内のSNPの一致度)
                        {
                            flagclust[i, j] = true;
                            flagclust[j, i] = true;
                        }
                    }
                }
            }
            List<int> maxmember = new List<int>();
            //最大クラスターのSNP群を見つける
            for (int i = 0; i < V_SNP数; i++)
            {
                if (!flagid[i])
                {
                    flagid[i] = true;
                    List<int> temp = new List<int>();
                    temp.Add(i);
                    //二次元のフラグ表の中で、iが入っている塊を抽出する
                    temp.AddRange(searchclust(flagclust, flagid, i));
                    if (temp.Count >= maxmember.Count)
                    {
                        maxmember = temp;
                    }
                }
            }
            int max = maxmember.Count;
            //Console.WriteLine(maxmember.Count);

            List<int[]> tempcdata = new List<int[]>();
            foreach (int i in maxmember)
            {
                tempcdata.Add(I_SNPごとのジェノタイプ一覧[i]);
            }
            int[] cons = get_consensus(tempcdata, rateOfNotNA);
            int tempa = maxmember.Min(x=>I_POS一覧[x]);
            int tempb = maxmember.Max(x=>I_POS一覧[x]);
            //最大クラスターについてのみブロックとして出力する
            writer.Write("sites: " + V_SNP数 + " in_cluster: " + maxmember.Count + " #"+maxmember.Min(x=>I_POS一覧[x]) +"#"+maxmember.Max(x=>I_POS一覧[x])
                             + "\t" + contig + "\t" + I_POS一覧[maxmember[0]]);
            foreach (int s in cons)
            {
                writer.Write("\t" + s);
            }
            writer.WriteLine("");
        }

        /// <summary>
        /// 二次元のフラグ表の中で、keyが入っている塊を抽出する
        /// </summary>
        /// <param name="flagclust">二次元のフラグ表</param>
        /// <param name="flagid">既に他の塊に入っているかどうか</param>
        /// <param name="key">検索するID</param>
        /// <returns></returns>
        public List<int> searchclust(bool[,] flagclust, bool[] flagid, int key)
        {
            List<int> result = new List<int>();
            for (int i = 0; i < flagid.Length; i++)
            {
                if (!flagid[i] && flagclust[key, i])
                {
                    flagid[i] = true;
                    result.Add(i);
                    result.AddRange(searchclust(flagclust, flagid, i));
                }
            }
            return result;
        }

        public int[] get_consensus(List<int[]> cdata, double rateOfNotNA)
        {
            if (cdata.Count == 0) { return null; }
            int[] baseseq = cdata[0];
            int[,] basecnt = new int[baseseq.Length, 2];
            calc_cons(baseseq, 1, basecnt);
            for (int i = 1; i < cdata.Count; i++)
            {
                distance temp = get_dist(baseseq, cdata[i], rateOfNotNA);
                calc_cons(cdata[i], temp.order, basecnt);
            }
            for (int i = 0; i < baseseq.Length; i++)
            {
                if (basecnt[i, 0] == basecnt[i, 1])
                {
                    baseseq[i] = -1;
                }
                else if (basecnt[i, 0] > basecnt[i, 1])
                {
                    baseseq[i] = 0;
                }
                else
                {
                    baseseq[i] = 1;
                }
            }
            return baseseq;
        }
        public void calc_cons(int[] seqs, int ord, int[,] basecnt)
        {
            for (int i = 0; i < seqs.Length; i++)
            {
                if ((ord == 1 && seqs[i] == 0) || (ord == -1 && seqs[i] == 1))
                {
                    basecnt[i, 0]++;
                }
                else if ((ord == 1 && seqs[i] == 1) || (ord == -1 && seqs[i] == 0))
                {
                    basecnt[i, 1]++;
                }
            }
        }
        public string arr2str(int[] x)
        {
            string str = x[0].ToString();
            for (int i = 1; i < x.Length; i++)
            {
                str = str + "\t" + x[i].ToString();
            }
            return str;
        }
    }
}