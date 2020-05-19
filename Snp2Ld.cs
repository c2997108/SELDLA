using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SELDLA
{
    class Snp2Ld
    {
        public Snp2Ld(){
            initialize();
        }

        public string oldchr;
        public int num;
        public List<int[]> data;
        public List<int> pos;
        public int th_r = 10000;
        public Encoding Enc = Encoding.GetEncoding("UTF-8");
        public StreamWriter writer;
        //public string inputvcf;
        //public string output;
        //public double opt_balance=0.1;
        public double opt_ld_match_rate = 0.9;
        public struct distance
        {
            public int order;
            public double match_rate;
        }
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
            oldchr = "";
            num = 0;
            data = new List<int[]>();
            pos = new List<int>();
        }
        public void run(string inputtxt, double opt_nc, int opt_r, double rateOfNotNA)
        {
            opt_ld_match_rate=opt_nc;
            th_r=opt_r;
            string output=inputtxt+".ld";
            int counter = 0;
            string line;

            // Read the file and display it line by line.  
            System.IO.StreamReader file = new System.IO.StreamReader(inputtxt);
            while ((line = file.ReadLine()) != null)
            {
                counter++;
                //System.Console.WriteLine(line);
                if (counter == 1)
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
                    int num_people = values.Length - 2;
                    if (oldchr != values[0] && oldchr != "" && num > 0)
                    {
                        ldsearch(oldchr, rateOfNotNA);
                        initialize();
                    }
                    oldchr = values[0];
                    num++;
                    int[] tempdata = new int[num_people];
                    for (int i = 2; i < values.Length; i++)
                    {
                        tempdata[i - 2] = Int32.Parse(values[i]);
                    }
                    data.Add(tempdata);
                    pos.Add(Int32.Parse(values[1]));
                }
            }
            if (num > 0)
            {
                ldsearch(oldchr, rateOfNotNA);
            }

            file.Close();
            writer.Close();
            System.Console.WriteLine("There were {0} lines.", counter);
            // Suspend the screen.  
            // System.Console.ReadLine();
        }

        public void ldsearch(string chr, double rateOfNotNA)
        {
            int cpos = 0;
            int cnum = 0;
            List<int[]> tempdata = new List<int[]>();
            List<int> temppos = new List<int>();
            for (int i = 0; i < num; i++)
            {
                if (cpos == 0) { cpos = pos[i]; }
                if (pos[i] <= cpos + th_r)
                {
                }
                else
                {
                    clustsearch(tempdata, temppos, chr, rateOfNotNA);
                    cpos = 0;
                    cnum = 0;
                    tempdata = new List<int[]>();
                    temppos = new List<int>();
                    cpos = pos[i];
                }
                cnum++;
                tempdata.Add(data[i]);
                temppos.Add(pos[i]);
            }
            clustsearch(tempdata, temppos, chr, rateOfNotNA);
        }

        public void clustsearch(List<int[]> cdata, List<int> cpos, string chr, double rateOfNotNA)
        {
            int numdata = cpos.Count;
            double[,] distmatrix = new double[numdata, numdata];
            int[,] corder = new int[numdata, numdata];
            bool[,] flagclust = new bool[numdata, numdata];
            bool[] flagid = new bool[numdata];
            for (int i = 0; i < numdata; i++)
            {
                for (int j = 0; j < numdata; j++)
                {
                    // i==jの場合は0で初期化されている
                    if (i > j)
                    {
                        distance tempdist = get_dist(cdata[i], cdata[j], rateOfNotNA);
                        distmatrix[i, j] = tempdist.match_rate;
                        distmatrix[j, i] = tempdist.match_rate;
                        corder[i, j] = tempdist.order;
                        corder[j, i] = tempdist.order;
                        if (tempdist.match_rate >= opt_ld_match_rate)
                        {
                            flagclust[i, j] = true;
                            flagclust[j, i] = true;
                        }
                    }
                }
            }
            List<int> maxmember = new List<int>();
            for (int i = 0; i < numdata; i++)
            {
                if (!flagid[i])
                {
                    flagid[i] = true;
                    List<int> temp = new List<int>();
                    temp.Add(i);
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
                tempcdata.Add(cdata[i]);
            }
            int[] cons = get_consensus(tempcdata, rateOfNotNA);
            //Console.Write("sites: " + numdata + " in_cluster: " + maxmember.Count + "\t" + chr + "\t" + cpos[maxmember[0]]);
            int tempa = maxmember.Min(x=>cpos[x]);
            int tempb = maxmember.Max(x=>cpos[x]);
            writer.Write("sites: " + numdata + " in_cluster: " + maxmember.Count + " #"+maxmember.Min(x=>cpos[x]) +"#"+maxmember.Max(x=>cpos[x])
                             + "\t" + chr + "\t" + cpos[maxmember[0]]);
            foreach (int s in cons)
            {
                //Console.Write("\t" + s);
                writer.Write("\t" + s);
            }
            //Console.WriteLine("");
            writer.WriteLine("");
        }
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