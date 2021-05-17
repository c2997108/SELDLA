using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SELDLA{
        public class SortableMatch : System.IComparable
    {
        public int V_contig1と2の長さ合計;
        public double V_一致率;
        public string contig1;
        public string contig2;
        public string start_end_1;
        public string start_end_2;

        public SortableMatch(int V_contig1と2の長さ合計, double V_一致率, string contig1, string contig2, string start_end_1, string start_end_2){
            this.V_contig1と2の長さ合計=V_contig1と2の長さ合計;
            this.V_一致率=V_一致率;
            this.contig1=contig1;
            this.contig2=contig2;
            this.start_end_1=start_end_1;
            this.start_end_2=start_end_2;
        }
        //自分自身がobjより小さいときはマイナスの数、大きいときはプラスの数、
        //同じときは0を返す
        public int CompareTo(object obj)
        {
            //nullより大きい
            if (obj == null)
            {
                return 1;
            }

            //違う型とは比較できない
            if (this.GetType() != obj.GetType())
            {
                throw new ArgumentException("別の型とは比較できません。", "obj");
            }
            //このクラスが継承されることが無い（構造体など）ならば、次のようにできる
            //if (!(other is TestClass)) { }

            if(this.V_一致率>((SortableMatch)obj).V_一致率){
                return 1;
            }else if(this.V_一致率<((SortableMatch)obj).V_一致率){
                return -1;
            }else{
                if(this.V_contig1と2の長さ合計>((SortableMatch)obj).V_contig1と2の長さ合計){
                    return 1;
                }else if(this.V_contig1と2の長さ合計<((SortableMatch)obj).V_contig1と2の長さ合計){
                    return -1;
                }else{
                    return 0;
                }
            }
        }
    }
}