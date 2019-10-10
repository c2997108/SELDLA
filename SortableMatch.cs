using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SELDLA{
        public class SortableMatch : System.IComparable
    {
        public int length;
        public double matchrate;
        public string chr1;
        public string chr2;
        public string pos1;
        public string pos2;

        public SortableMatch(int length, double matchrate, string chr1, string chr2, string pos1, string pos2){
            this.length=length;
            this.matchrate=matchrate;
            this.chr1=chr1;
            this.chr2=chr2;
            this.pos1=pos1;
            this.pos2=pos2;
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

            if(this.matchrate>((SortableMatch)obj).matchrate){
                return 1;
            }else if(this.matchrate<((SortableMatch)obj).matchrate){
                return -1;
            }else{
                if(this.length>((SortableMatch)obj).length){
                    return 1;
                }else if(this.length<((SortableMatch)obj).length){
                    return -1;
                }else{
                    return 0;
                }
            }
        }
    }
}