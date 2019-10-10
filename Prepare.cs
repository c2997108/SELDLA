using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Mono.Options;
using Bio;
using Bio.IO;

namespace SELDLA
{
    class Prepare
    {
        public Prepare()
        {
            initialize();
        }
        public void initialize()
        {
        }
        public void cleanupVcf(string filename, int opt_dp, int opt_gq, double opt_nonzerorate, string opt_o)
        {
            StreamWriter writer = new StreamWriter(opt_o + "_clean.txt");
            // System.IO.StreamReader file = new System.IO.StreamReader(filename);
            // int counter = 0;
            // string line;
            // //List<string> lines=new List<string>();
            // int outline=0;
            // while ((line = file.ReadLine()) != null)
            // {
            //     if(!line.StartsWith("##")){
            //         counter++;
            //         if(counter==1){
            //             //Console.WriteLine(line);
            //             writer.WriteLine(line);
            //         }else{
            //             //lines.Add(line);
            //             string str = cleanupRow(line,opt_dp, opt_gq, opt_nonzerorate);
            //             if(str!=null){
            //                 outline++;
            //                 writer.WriteLine(str);
            //             }
            //         }
            //     }
            // }
            // file.Close();
            // Console.WriteLine("cleaned datas: "+outline);

            // ParallelQuery<string> cleaned=lines
            //                                 .AsParallel()
            //                                 .AsOrdered()
            //                                 .Select(f => cleanupRow(f, opt_dp, opt_gq, opt_nonzerorate));
            List<string> cleaned = GetRecordsFromFile(filename).AsParallel().AsOrdered().Select(x => cleanupRow(x, opt_dp, opt_gq, opt_nonzerorate)).ToList();
            Console.WriteLine("input datas: " + cleaned.Count());
            cleaned = cleaned.Where(x => x != null).OrderBy(x => x.Split("\t")[0]).ThenBy(x => ParseInt(x.Split("\t")[1])).ToList();
            Console.WriteLine("cleaned datas: " + cleaned.Count());
            foreach (string str in cleaned)
            {
                //if (str != null)
                //{
                //Console.WriteLine(str);
                writer.WriteLine(str);
                //}
            }
            writer.Close();
        }

        private static int ParseInt(string str)
        {
            int result;

            int.TryParse(str, out result);

            return result;
        }

        static IEnumerable<string> GetRecordsFromFile(string filename)
        {
            using (var streamReader = new StreamReader(filename))
            {
                string str;
                while ((str = streamReader.ReadLine()) != null)
                {
                    yield return str;
                }
            }
        }

        public static string cleanupRow(string line, int opt_dp, int opt_gq, double opt_nonzerorate)
        {
            bool viewflag = true;
            string result = "";
            if (line.StartsWith("##"))
            {
                viewflag = false;
            }
            else if (line.StartsWith("#"))
            {
                result = line;
            }
            else
            {
                string[] arr1 = line.Split("\t");
                result = arr1[0];
                if (arr1[3].Length == 1 && arr1[4].Length == 1)
                {
                    string[] descript = arr1[8].Split(":");
                    int pos_dp = 0;
                    int pos_gq = 0;
                    for (int i = 1; i < 7; i++)
                    {
                        result += "\t" + arr1[i];
                    }
                    result += "\t.\t" + arr1[8];
                    for (int i = 0; i < descript.Length; i++)
                    {
                        if (descript[i] == "DP")
                        {
                            pos_dp = i;
                            //Console.WriteLine("DP" + i);
                        }
                        else if (descript[i] == "GQ")
                        {
                            pos_gq = i;
                            //Console.WriteLine("GQ" + i);
                        }
                    }
                    int numofnonzero = 0;
                    for (int i = 9; i < arr1.Length; i++)
                    {
                        string[] value_pearson = arr1[i].Split(":");
                        if ((pos_dp == 0 || (value_pearson.Length > pos_dp && value_pearson[pos_dp] != "." && Int32.Parse(value_pearson[pos_dp]) >= opt_dp)) &&
                            (pos_gq == 0 || (value_pearson.Length > pos_gq && value_pearson[pos_gq] != "." && Int32.Parse(value_pearson[pos_gq]) >= opt_gq)))
                        {
                            string tempres = "";
                            if (value_pearson[0] == "0/0")
                            {
                                tempres = "0";
                                numofnonzero++;
                            }
                            else if (value_pearson[0] == "0/1" || value_pearson[0] == "0|1" || value_pearson[0] == "1|0")
                            {
                                tempres = "1";
                                numofnonzero++;
                            }
                            else if (value_pearson[0] == "1/1" || value_pearson[0] == "1|1")
                            {
                                tempres = "2";
                                numofnonzero++;
                            }
                            else
                            {
                                tempres = "-1";
                            }
                            result += "\t" + tempres;
                        }
                        else
                        {
                            result += "\t-1";
                        }
                    }
                    if (numofnonzero < (arr1.Length - 9) * opt_nonzerorate)
                    {
                        viewflag = false;
                    }
                }
                else
                {
                    viewflag = false;
                }

            }
            if (viewflag)
            {
                return result;
            }
            else
            {
                return null;
            }
        }

        public void splitVcf(string input_clean, string opt_o, string inputfamily, double opt_p, double opt_b, string mode)
        {
            StreamReader file = new StreamReader(inputfamily);
            string line;
            List<string[]> families = new List<string[]>();
            while ((line = file.ReadLine()) != null)
            {
                string[] temp = line.Split("\t");
                if ((mode == "crossbreed" || mode == "duploid") && temp.Length >= 3)
                {
                    families.Add(temp);
                }
                else if (mode == "haploid" && temp.Length >= 2)
                {
                    families.Add(temp);
                }
            }
            file.Close();
            Console.WriteLine("families: " + families.Count);

            StreamReader vcffile = new StreamReader(input_clean);
            List<string[]> lines = new List<string[]>();
            string[] header = null;
            int counter = 0;
            while ((line = vcffile.ReadLine()) != null)
            {
                counter++;
                if (counter == 1)
                {
                    header = line.Split("\t");
                }
                else
                {
                    lines.Add(line.Split("\t"));
                }
            }
            vcffile.Close();
            Console.WriteLine("cleaned datas: " + lines.Count);
            lines = lines.OrderBy(x => x[0]).ThenBy(x => Int32.Parse(x[1])).ToList();

            int num_fam = 0;
            foreach (string[] fam in families)
            {
                num_fam++;
                Dictionary<int, int> idord = new Dictionary<int, int>();
                for (int j = 0; j < fam.Length; j++)
                {
                    if ((mode == "crossbreed" || mode == "duploid") && (j == 0 || j == 1) && fam[j] == "-1")
                    {
                        idord.Add(j, -1);
                    }
                    else if (mode == "haploid" && j == 0 && fam[j] == "-1")
                    {
                        idord.Add(j, -1);
                    }
                    else
                    {
                        for (int i = 9; i < header.Length; i++)
                        {
                            if (fam[j] == header[i])
                            {
                                idord.Add(j, i);
                            }
                        }
                        if (!idord.ContainsKey(j))
                        {
                            throw new System.Exception("The corresponding ID of " + fam[j] + " does not exist in the VCF.");
                        }
                    }
                }
                StreamWriter writer = new StreamWriter(opt_o + "_split_" + num_fam + ".txt");
                writer.Write("#chr\tpos");
                if (mode == "crossbreed" || mode == "duploid")
                {
                    for (int j = 2; j < fam.Length; j++)
                    {
                        writer.Write("\t" + header[idord[j]]);
                    }
                }
                else
                { //haploid
                    for (int j = 1; j < fam.Length; j++)
                    {
                        writer.Write("\t" + header[idord[j]]);
                    }
                }
                writer.WriteLine("");

                ParallelQuery<string> datas;
                if (mode == "crossbreed")
                {
                    datas = lines
                                .AsParallel()
                                .AsOrdered()
                                .Select(f => splitVcfRowCrossHap(f, idord, opt_p, opt_b));
                }
                else if (mode == "duploid")
                {
                    datas = lines
                                .AsParallel()
                                .AsOrdered()
                                .Select(f => splitVcfRowDup(f, idord, opt_p, opt_b));
                }
                else
                {
                    datas = lines
                                .AsParallel()
                                .AsOrdered()
                                .Select(f => splitVcfRowHap(f, idord, opt_p, opt_b));
                }
                int outline = 0;
                foreach (string str in datas)
                {
                    if (str != null)
                    {
                        outline++;
                        //Console.WriteLine(str);
                        writer.WriteLine(str);
                    }
                }
                writer.Close();
                Console.WriteLine("family " + num_fam + " datas: " + outline);
            }

        }

        public static string splitVcfRowCrossHap(string[] vals, Dictionary<int, int> idord, double opt_p, double opt_b)
        {
            string result = "";
            bool viewflag = true;
            if ((idord[0]!=-1 && vals[idord[0]] != "1") || (idord[1]!=-1 && vals[idord[1]] != "-1")) { viewflag = false; }
            if (viewflag)
            {
                int n0 = 0;
                int n1 = 0;
                int n2 = 0;
                int n = 0;
                for (int i = 2; i < idord.Count; i++)
                {
                    n++;
                    if (vals[idord[i]] == "0")
                    {
                        n0++;
                    }
                    else if (vals[idord[i]] == "1")
                    {
                        n1++;
                    }
                    else if (vals[idord[i]] == "2")
                    {
                        n2++;
                    }
                }
                if ((n0 + n2) < n * opt_p || n0 < n * opt_b || n2 < n * opt_b)
                {
                    viewflag = false;
                }
                else
                {
                    result = vals[0] + "\t" + vals[1];
                    for (int i = 2; i < idord.Count; i++)
                    {
                        if (vals[idord[i]] == "0" || vals[idord[i]] == "-1")
                        {
                            result += "\t" + vals[idord[i]];
                        }
                        else if (vals[idord[i]] == "2")
                        {
                            result += "\t1";
                        }
                        else
                        {
                            result += "\t-1";
                        }
                    }
                }
            }
            if (!viewflag)
            {
                result = null;
            }
            return result;
        }

        public static string splitVcfRowHap(string[] vals, Dictionary<int, int> idord, double opt_p, double opt_b)
        {
            string result = "";
            bool viewflag = true;
            if (idord[0]!=-1 && vals[idord[0]] != "1") { viewflag = false; }
            if (viewflag)
            {
                int n0 = 0;
                int n1 = 0;
                int n2 = 0;
                int n = 0;
                for (int i = 1; i < idord.Count; i++)
                {
                    n++;
                    if (vals[idord[i]] == "0")
                    {
                        n0++;
                    }
                    else if (vals[idord[i]] == "1")
                    {
                        n1++;
                    }
                    else if (vals[idord[i]] == "2")
                    {
                        n2++;
                    }
                }
                if ((n0 + n2) < n * opt_p || n0 < n * opt_b || n2 < n * opt_b)
                {
                    viewflag = false;
                }
                else
                {
                    result = vals[0] + "\t" + vals[1];
                    for (int i = 1; i < idord.Count; i++)
                    {
                        if (vals[idord[i]] == "0" || vals[idord[i]] == "-1")
                        {
                            result += "\t" + vals[idord[i]];
                        }
                        else if (vals[idord[i]] == "2")
                        {
                            result += "\t1";
                        }
                        else
                        {
                            result += "\t-1";
                        }
                    }
                }
            }
            if (!viewflag)
            {
                result = null;
            }
            return result;
        }

        public static string splitVcfRowDup(string[] vals, Dictionary<int, int> idord, double opt_p, double opt_b)
        {
            string result = "";
            bool viewflag = true;
            if ((idord[0]!=-1 && vals[idord[0]] != "1") || (idord[1]!=-1 && !(vals[idord[1]] == "0" || vals[idord[1]] == "2"))) { viewflag = false; }
            if (viewflag)
            {
                int n0 = 0;
                int n1 = 0;
                int n2 = 0;
                int n = 0;
                for (int i = 2; i < idord.Count; i++)
                {
                    n++;
                    if (vals[idord[i]] == "0")
                    {
                        n0++;
                    }
                    else if (vals[idord[i]] == "1")
                    {
                        n1++;
                    }
                    else if (vals[idord[i]] == "2")
                    {
                        n2++;
                    }
                }
                if (n1 < n * opt_b
                 || (vals[idord[1]] == "2" && (n2 < n * opt_b || (n1 + n2) < n * opt_p))
                 || (vals[idord[1]] == "0" && (n0 < n * opt_b || (n1 + n0) < n * opt_p)))
                {
                    viewflag = false;
                }
                else
                {
                    result = vals[0] + "\t" + vals[1];
                    for (int i = 2; i < idord.Count; i++)
                    {
                        if (vals[idord[i]] == "1" || vals[idord[i]] == "-1")
                        {
                            result += "\t" + vals[idord[i]];
                        }
                        else if ((vals[idord[1]] == "2" && vals[idord[i]] == "2") || (vals[idord[1]] == "0" && vals[idord[i]] == "0"))
                        {
                            result += "\t0";
                        }
                        else
                        {
                            result += "\t-1";
                        }
                    }
                }
            }
            if (!viewflag)
            {
                result = null;
            }
            return result;
        }
    }
}
