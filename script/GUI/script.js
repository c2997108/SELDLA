
var chr2nd = [];
var chrorig = [];
var chrorient = [];
var dataphase = [];
var distphase = [];
var markerpos = [];
// GPU is a constructor and namespace for browser
const gpu = new GPU();
//GPUを使ってマーカー間の一致率計算
function calcdist(a) {
  const multiplyMatrix = gpu.createKernel(function (a, size) {
    let sum = 0;
    let sum2 = 0;
    let n = 0;
    for (let i = 0; i < size; i++) {
      if (a[this.thread.y][i] != 0 && a[this.thread.x][i] != 0) {
        n++;
        //if(a[this.thread.y][i] * b[i][this.thread.x]>0.5){sum++}
        //if(a[this.thread.y][i] * b[i][this.thread.x]<-0.5){sum2++}
        if (a[this.thread.y][i] == a[this.thread.x][i]) { sum++ }
        if (a[this.thread.y][i] == -a[this.thread.x][i]) { sum2++ }
      }
    }
    if (sum > sum2) { return (sum + 0.0001) / n } else { return (sum2 + 0.0001) / n } //何故かNVIDIA RTX2060ではちょうど0.5になる場合0.9999999403953552として返ってくる
    //if(sum > sum2){return n}else{return n}
  }, { loopMaxIterations: a.length }).setOutput([a.length, a.length]);

  //return multiplyMatrix(a, transpose(a), a[0].length)
  return multiplyMatrix(a, a[0].length);
}

//行列を転置させる関数定義
//const transpose = a => a[0].map((_, c) => a.map(r => r[c]));
function transpose(a) {
  let resarray = [];
  for (let i = 0; i < a[0].length; i++) {
    let temparray = [];
    for (let j = 0; j < a.length; j++) {
      temparray.push(a[j][i]);
    }
    resarray.push(temparray);
  }
  return resarray;
}



var lines;
window.onload = function () {
  var file = document.querySelector('#getfile');
  file.onchange = function () {
    const spinner = document.getElementById('loading');
    spinner.classList.remove('loaded');
    setTimeout(function () {
      var fileList = file.files;
      var reader = new FileReader();
      reader.readAsText(fileList[0]);
      //読み込み後表示
      reader.onload = function () {
        chr2nd = [];
        chrorig = [];
        chrorient = [];
        dataphase = [];
        distphase = [];
        markerpos = [];
        lines = this.result.split(/\r\n|\n/);
        g.clear();
        g.lineStyle(0, 0xFFFF00);
        for (let j = 0; j < lines.length; j++) {
          if (lines[j] !== "") { //最後改行してある普通のファイルなら空の値が入ってくる
            //console.log(j+": "+lines[j]);
            if (lines[j].substr(0, 1) !== "#") {
              let items = lines[j].split(/\t/);
              chr2nd.push(items[0]);
              chrorig.push(items[4]);
              if ((items[2] === "+" && items[3] === "+") || (items[2] === "-" && items[3] === "-")) {
                chrorient.push("+");
              } else if ((items[2] === "+" && items[3] === "-") || (items[2] === "-" && items[3] === "+")) {
                chrorient.push("-");
              } else { chrorient.push("na"); }
              let tempdata = [];
              for (let i = 6; i < items.length; i++) { if (items[i] === "1") { tempdata.push(1) } else if (items[i] === "0") { tempdata.push(-1) } else { tempdata.push(0) } }
              dataphase.push(tempdata);
              markerpos.push(items[5]);
            }
          }
        }
        distphase = calcdist(dataphase);
        //console.log(distphase);
        //console.log(dataphase);
        for (let j = 0; j < distphase.length; j++) {
          for (let i = 0; i < distphase.length; i++) {
            g.beginFill(16777215 - (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)) - 256 * (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)));
            g.drawRect(j, i, 1, 1);
            g.endFill();
          }
        }
        savedata();
        spinner.classList.add('loaded');
      }
    }, 1000);
  };



  fetch('http://suikou.fs.a.u-tokyo.ac.jp/test/seldla/savedate_7.txt')
    .then(response => response.text())
    .then(data => {
      lines = data.split(/\r\n|\n/);
      g.clear();
      g2.clear();
      g.lineStyle(0, 0xFFFF00);
      for (let j = 0; j < lines.length; j++) {
        if (lines[j] !== "") { //最後改行してある普通のファイルなら空の値が入ってくる
          //console.log(j+": "+lines[j]);
          if (lines[j].substr(0, 1) !== "#") {
            let items = lines[j].split(/\t/);
            chr2nd.push(items[0]);
            chrorig.push(items[4]);
            if ((items[2] === "+" && items[3] === "+") || (items[2] === "-" && items[3] === "-")) {
              chrorient.push("+");
            } else if ((items[2] === "+" && items[3] === "-") || (items[2] === "-" && items[3] === "+")) {
              chrorient.push("-");
            } else { chrorient.push("na"); }
            let tempdata = [];
            for (let i = 6; i < items.length; i++) { if (items[i] === "1") { tempdata.push(1) } else if (items[i] === "0") { tempdata.push(-1) } else { tempdata.push(0) } }
            dataphase.push(tempdata);
            markerpos.push(items[5]);
          }
        }
      }
      //console.log(dataphase);
      distphase = calcdist(dataphase);
      for (let j = 0; j < distphase.length; j++) {
        for (let i = 0; i < distphase.length; i++) {
          g.beginFill(16777215 - (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)) - 256 * (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)));
          g.drawRect(j, i, 1, 1);
          g.endFill();
        }
      }
      savedata();
      const spinner = document.getElementById('loading');
      spinner.classList.add('loaded');

    });

}










let app = new PIXI.Application({ width: window.innerWidth - 30, height: window.innerHeight - 80 });
app.stage.scale.x = 1.1;
app.stage.scale.y = 1.1;
let pixijscanvas = document.body.appendChild(app.view);
var g = new PIXI.Graphics();
var g2 = new PIXI.Graphics();

var markN = 0;
var markM = 0;
var markNstart = [0, 0, 0];
var markNend = [0, 0, 0];
var markMstart1 = [0, 0, 0];
var markMend1 = [0, 0, 0];
var markMstart2 = [0, 0, 0];
var markMend2 = [0, 0, 0];

g.lineStyle(0, 0xFFFF00);

app.stage.addChild(g);

//const message = new PIXI.Text("Loading...");
//message.x = 200;
//message.y = 200;
//message.style.fill = "white";
//app.stage.addChild(message);
//message.destroy(true);

//g2.lineStyle(1, 0xFFFF00);
//g2.drawRect(0,0,100,100);
app.stage.addChild(g2);

function zoom(s, x, y) {

  s = s > 0 ? 1.1 : 0.9;
  var worldPos = { x: (x - app.stage.x) / app.stage.scale.x, y: (y - app.stage.y) / app.stage.scale.y };
  var newScale = { x: app.stage.scale.x * s, y: app.stage.scale.y * s };

  var newScreenPos = { x: (worldPos.x) * newScale.x + app.stage.x, y: (worldPos.y) * newScale.y + app.stage.y };

  app.stage.x -= (newScreenPos.x - x);
  app.stage.y -= (newScreenPos.y - y);
  app.stage.scale.x = newScale.x;
  app.stage.scale.y = newScale.y;
};

var lastPos = null
$(pixijscanvas)
  .mousewheel(function (e) {
    zoom(e.deltaY, e.offsetX, e.offsetY)
  }).mousedown(function (e) {
    lastPos = { x: e.offsetX, y: e.offsetY };
  }).mouseup(function (event) {
    lastPos = null;
  }).mousemove(function (e) {
    if (lastPos) {
      app.stage.x += (e.offsetX - lastPos.x);
      app.stage.y += (e.offsetY - lastPos.y);
      lastPos = { x: e.offsetX, y: e.offsetY };
    }

    let tempX1 = parseInt((e.offsetX - app.stage.x) / app.stage.scale.x);
    let tempY1 = parseInt((e.offsetY - app.stage.y) / app.stage.scale.y);
    if (tempX1 < tempY1) { let tempX1temp = tempX1; tempX1 = tempY1; tempY1 = tempX1temp }
    $('#x1').text(tempX1);
    $('#y1').text(tempY1);
    if (tempX1 >= 0 && tempX1 <= chr2nd.length && tempY1 >= 0 && tempY1 <= chr2nd.length) {
      $('#mytooltip1').html("[Up side]<br>"
        + "Block number: " + tempY1 + "<br>"
        + "Chr: " + chr2nd[tempY1] + "<br>"
        + "Contig name: " + chrorig[tempY1] + "<br>"
        + "Contig orientation: " + chrorient[tempY1] + "<br>"
        + "Marker position: " + markerpos[tempY1] + "<br>"
        + "[Down side]<br>"
        + "Block number: " + tempX1 + "<br>"
        + "Chr: " + chr2nd[tempX1] + "<br>"
        + "Contig name: " + chrorig[tempX1] + "<br>"
        + "Contig orientation: " + chrorient[tempX1] + "<br>"
        + "Marker position: " + markerpos[tempX1] + "<br>"
        + "<br>Phase identity: " + distphase[tempX1][tempY1]);
      g2.clear();
      let tempReplace1 = chr2nd[tempX1];
      let tempReplace2 = chr2nd[tempY1];
      let tempReplace3 = chrorig[tempX1];
      let tempReplace4 = chrorig[tempY1];
      let tempStart = [-1, -1, -1, -1];
      let tempSize = [0, 0, 0, 0];
      for (let i = 0; i < chr2nd.length; i++) {
        if (chr2nd[i] === tempReplace1) {
          if (tempStart[0] === -1) { tempStart[0] = i };
          tempSize[0]++;
        }
        if (chr2nd[i] === tempReplace2) {
          if (tempStart[1] === -1) { tempStart[1] = i };
          tempSize[1]++;
        }
        if (chrorig[i] === tempReplace3) {
          if (tempStart[2] === -1) { tempStart[2] = i };
          tempSize[2]++;
        }
        if (chrorig[i] === tempReplace4) {
          if (tempStart[3] === -1) { tempStart[3] = i };
          tempSize[3]++;
        }
      }
      //g2.drawRect(parseInt((e.offsetX-app.stage.x)/app.stage.scale.x),parseInt((e.offsetY-app.stage.y)/app.stage.scale.y),1,1);
      g2.lineStyle(0.7, 0xFFFF00);
      g2.drawRect(tempStart[0], tempStart[0], tempSize[0], tempSize[0]);
      g2.drawRect(tempStart[1], tempStart[1], tempSize[1], tempSize[1]);
      g2.lineStyle(0.7, 0x00FFFF);
      g2.drawRect(tempStart[2], tempStart[2], tempSize[2], tempSize[2]);
      g2.drawRect(tempStart[3], tempStart[3], tempSize[3], tempSize[3]);

      if (markN === 1) {
        let tempNstart;
        let tempNend;
        if (tempStart[3] < markNstart[0]) {
          tempNstart = tempStart[3];
          tempNend = markNend[0];
        } else {
          tempNstart = markNstart[0];
          tempNend = tempStart[3] + tempSize[3] - 1;
        }
        //console.log(tempNstart,tempNstart,tempNend-tempNstart+1,tempNend-tempNstart+1);
        g2.lineStyle(0.7, 0xFF0000);
        g2.drawRect(tempNstart, tempNstart, tempNend - tempNstart + 1, tempNend - tempNstart + 1);
      }

      if (markM === 1) {
        let tempMstart1;
        let tempMend1;
        if (tempStart[3] < markMstart1[0]) {
          tempMstart1 = tempStart[3];
          tempMend1 = markMend1[0];
        } else {
          tempMstart1 = markMstart1[0];
          tempMend1 = tempStart[3] + tempSize[3] - 1;
        }
        g2.lineStyle(0.7, 0xFF0000);
        g2.drawRect(tempMstart1, tempMstart1, tempMend1 - tempMstart1 + 1, tempMend1 - tempMstart1 + 1);
      } else if (markM === 2) {
        let tempMstart1 = markMstart1[2];
        let tempMend1 = markMend1[2];
        g2.lineStyle(0.7, 0xFF0000);
        g2.drawRect(tempMstart1, tempMstart1, tempMend1 - tempMstart1 + 1, tempMend1 - tempMstart1 + 1);
      } else if (markM === 3) {
        let tempMstart1 = markMstart1[2];
        let tempMend1 = markMend1[2];
        let tempMstart2;
        let tempMend2;
        if (tempStart[3] < markMstart2[0]) {
          tempMstart2 = tempStart[3];
          tempMend2 = markMend2[0];
        } else {
          tempMstart2 = markMstart2[0];
          tempMend2 = tempStart[3] + tempSize[3] - 1;
        }
        g2.lineStyle(0.7, 0xFF0000);
        g2.drawRect(tempMstart1, tempMstart1, tempMend1 - tempMstart1 + 1, tempMend1 - tempMstart1 + 1);
        g2.drawRect(tempMstart2, tempMstart2, tempMend2 - tempMstart2 + 1, tempMend2 - tempMstart2 + 1);
      }

    }
  });

$(window).keydown(function (e) {
  console.log(e.keyCode);

  let tempX1 = parseInt($('#x1').text());
  let tempY1 = parseInt($('#y1').text());
  if (tempX1 < tempY1) { let tempX1temp = tempX1; tempX1 = tempY1; tempY1 = tempX1temp }
  let tempReplace1 = chr2nd[tempX1];
  let tempReplace2 = chr2nd[tempY1];
  let tempReplace3 = chrorig[tempX1];
  let tempReplace4 = chrorig[tempY1];
  if (e.keyCode === 82) {
    console.log("r");
    const spinner = document.getElementById('loading');
    spinner.classList.remove('loaded');
    setTimeout(function () {
      let tempchr2nd = [];
      let tempchrorig = [];
      let tempchrorient = [];
      let tempdataphase = [];
      let tempmarkerpos = [];
      let flag1 = 0;
      let flag2 = 0;
      for (let i = 0; i < chr2nd.length; i++) {
        if (chr2nd[i] !== tempReplace1) {
          tempchr2nd.push(chr2nd[i]);
          tempchrorig.push(chrorig[i]);
          tempchrorient.push(chrorient[i]);
          tempdataphase.push(dataphase[i]);
          tempmarkerpos.push(markerpos[i]);
        } else if (chr2nd[i] === tempReplace1 && flag1 === 0) {
          flag1 = 1;
          for (let j = chr2nd.length - 1; j >= 0; j--) {
            if (chr2nd[j] === tempReplace1) {
              tempchr2nd.push(chr2nd[j]);
              tempchrorig.push(chrorig[j]);
              if (chrorient[j] === "+") {
                tempchrorient.push("-");
              } else if (chrorient[j] === "-") {
                tempchrorient.push("+");
              } else {
                tempchrorient.push("na");
              }
              tempdataphase.push(dataphase[j]);
              tempmarkerpos.push(markerpos[j]);
            }
          }
        }
      }
      chr2nd = tempchr2nd;
      chrorig = tempchrorig;
      chrorient = tempchrorient;
      dataphase = tempdataphase;
      //console.log(dataphase);
      distphase = calcdist(dataphase);
      markerpos = tempmarkerpos;
      savedata();

      g.clear();
      for (let j = 0; j < distphase.length; j++) {
        for (let i = 0; i < distphase.length; i++) {
          g.beginFill(16777215 - (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)) - 256 * (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)));
          g.drawRect(j, i, 1, 1);
          g.endFill();
        }
      }

      spinner.classList.add('loaded');
    }, 1000);
  }

  if (e.keyCode === 84) {
    console.log("t");
    const spinner = document.getElementById('loading');
    spinner.classList.remove('loaded');
    setTimeout(function () {
      if (tempReplace1 !== tempReplace2) {
        let tempchr2nd = [];
        let tempchrorig = [];
        let tempchrorient = [];
        let tempdataphase = [];
        let tempmarkerpos = [];
        let flag1 = 0;
        let flag2 = 0;
        for (let i = 0; i < chr2nd.length; i++) {
          if (chr2nd[i] !== tempReplace2 && chr2nd[i] !== tempReplace1) {
            tempchr2nd.push(chr2nd[i]);
            tempchrorig.push(chrorig[i]);
            tempchrorient.push(chrorient[i]);
            tempdataphase.push(dataphase[i]);
            tempmarkerpos.push(markerpos[i]);
          } else if (chr2nd[i] === tempReplace2 && flag2 === 0) {
            flag2 = 1;
            for (let j = 0; j < chr2nd.length; j++) {
              if (chr2nd[j] === tempReplace1) {
                tempchr2nd.push(chr2nd[j]);
                tempchrorig.push(chrorig[j]);
                tempchrorient.push(chrorient[j]);
                tempdataphase.push(dataphase[j]);
                tempmarkerpos.push(markerpos[j]);
              }
            }
          } else if (chr2nd[i] === tempReplace1 && flag1 === 0) {
            flag1 = 1;
            for (let j = 0; j < chr2nd.length; j++) {
              if (chr2nd[j] === tempReplace2) {
                tempchr2nd.push(chr2nd[j]);
                tempchrorig.push(chrorig[j]);
                tempchrorient.push(chrorient[j]);
                tempdataphase.push(dataphase[j]);
                tempmarkerpos.push(markerpos[j]);
              }
            }
          }
        }
        chr2nd = tempchr2nd;
        chrorig = tempchrorig;
        chrorient = tempchrorient;
        dataphase = tempdataphase;
        //console.log(dataphase);
        distphase = calcdist(dataphase);
        markerpos = tempmarkerpos;
        savedata();
      }

      g.clear();
      for (let j = 0; j < distphase.length; j++) {
        for (let i = 0; i < distphase.length; i++) {
          g.beginFill(16777215 - (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)) - 256 * (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)));
          g.drawRect(j, i, 1, 1);
          g.endFill();
        }
      }
      spinner.classList.add('loaded');
    }, 1000)
  }

  if (e.keyCode === 89) {
    console.log("y");
    const spinner = document.getElementById('loading');
    spinner.classList.remove('loaded');
    setTimeout(function () {
      let tempchr2nd = [];
      let tempchrorig = [];
      let tempchrorient = [];
      let tempdataphase = [];
      let tempmarkerpos = [];
      let flag1 = 0;
      let flag2 = 0;
      for (let i = 0; i < chr2nd.length; i++) {
        if (chrorig[i] !== tempReplace3) {
          tempchr2nd.push(chr2nd[i]);
          tempchrorig.push(chrorig[i]);
          tempchrorient.push(chrorient[i]);
          tempdataphase.push(dataphase[i]);
          tempmarkerpos.push(markerpos[i]);
        } else if (chrorig[i] === tempReplace3 && flag1 === 0) {
          flag1 = 1;
          for (let j = chr2nd.length - 1; j >= 0; j--) {
            if (chrorig[j] === tempReplace3) {
              tempchr2nd.push(chr2nd[j]);
              tempchrorig.push(chrorig[j]);
              if (chrorient[j] === "+") {
                tempchrorient.push("-");
              } else if (chrorient[j] === "-") {
                tempchrorient.push("+");
              } else {
                tempchrorient.push("na");
              }
              tempdataphase.push(dataphase[j]);
              tempmarkerpos.push(markerpos[j]);
            }
          }
        }
      }
      chr2nd = tempchr2nd;
      chrorig = tempchrorig;
      chrorient = tempchrorient;
      dataphase = tempdataphase;
      //console.log(dataphase);
      distphase = calcdist(dataphase);
      markerpos = tempmarkerpos;
      savedata();

      g.clear();
      for (let j = 0; j < distphase.length; j++) {
        for (let i = 0; i < distphase.length; i++) {
          g.beginFill(16777215 - (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)) - 256 * (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)));
          g.drawRect(j, i, 1, 1);
          g.endFill();
        }
      }
      spinner.classList.add('loaded');
    }, 1000)
  }

  if (e.keyCode === 85) {
    console.log("u");
    const spinner = document.getElementById('loading');
    spinner.classList.remove('loaded');
    setTimeout(function () {
      if (tempReplace3 !== tempReplace4) {
        let tempchr2nd = [];
        let tempchrorig = [];
        let tempchrorient = [];
        let tempdataphase = [];
        let tempmarkerpos = [];
        let flag1 = 0;
        let flag2 = 0;
        for (let i = 0; i < chr2nd.length; i++) {
          if (chrorig[i] !== tempReplace4 && chrorig[i] !== tempReplace3) {
            tempchr2nd.push(chr2nd[i]);
            tempchrorig.push(chrorig[i]);
            tempchrorient.push(chrorient[i]);
            tempdataphase.push(dataphase[i]);
            tempmarkerpos.push(markerpos[i]);
          } else if (chrorig[i] === tempReplace4 && flag2 === 0) {
            flag2 = 1;
            for (let j = 0; j < chr2nd.length; j++) {
              if (chrorig[j] === tempReplace3) {
                tempchr2nd.push(chr2nd[j]);
                tempchrorig.push(chrorig[j]);
                tempchrorient.push(chrorient[j]);
                tempdataphase.push(dataphase[j]);
                tempmarkerpos.push(markerpos[j]);
              }
            }
          } else if (chrorig[i] === tempReplace3 && flag1 === 0) {
            flag1 = 1;
            for (let j = 0; j < chr2nd.length; j++) {
              if (chrorig[j] === tempReplace4) {
                tempchr2nd.push(chr2nd[j]);
                tempchrorig.push(chrorig[j]);
                tempchrorient.push(chrorient[j]);
                tempdataphase.push(dataphase[j]);
                tempmarkerpos.push(markerpos[j]);
              }
            }
          }
        }
        chr2nd = tempchr2nd;
        chrorig = tempchrorig;
        chrorient = tempchrorient;
        dataphase = tempdataphase;
        //console.log(dataphase);
        distphase = calcdist(dataphase);
        markerpos = tempmarkerpos;
        savedata();
      }

      g.clear();
      for (let j = 0; j < distphase.length; j++) {
        for (let i = 0; i < distphase.length; i++) {
          g.beginFill(16777215 - (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)) - 256 * (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)));
          g.drawRect(j, i, 1, 1);
          g.endFill();
        }
      }
      spinner.classList.add('loaded');
    })
  }

  if (e.keyCode === 78) {
    console.log("n");
    if (markN === 0) {
      markN = 1;
      markNstart[0] = -1;
      let tempReplace3 = chrorig[tempX1];
      for (let i = 0; i < chrorig.length; i++) {
        if (chrorig[i] === tempReplace3) {
          if (markNstart[0] === -1) { markNstart[0] = i };
          markNend[0] = i;
        }
      }
    } else if (markN === 1) {
      const spinner = document.getElementById('loading');
      spinner.classList.remove('loaded');
      setTimeout(function () {
        markN = 0;
        markNstart[1] = -1;
        let tempReplace3 = chrorig[tempX1];
        for (let i = 0; i < chrorig.length; i++) {
          if (chrorig[i] === tempReplace3) {
            if (markNstart[1] === -1) { markNstart[1] = i };
            markNend[1] = i;
          }
        }
        if (markNstart[0] < markNstart[1]) {
          markNstart[2] = markNstart[0];
          markNend[2] = markNend[1];
        } else {
          markNstart[2] = markNstart[1];
          markNend[2] = markNend[0];
        }
        console.log(markNstart + ", " + markNend);
        g2.lineStyle(1, 0xFF0000);
        g2.drawRect(markNstart, markNstart, markNend - markNstart + 1, markNend - markNstart + 1);

        let tempchr2nd = [];
        let tempchrorig = [];
        let tempchrorient = [];
        let tempdataphase = [];
        let tempmarkerpos = [];
        let flag1 = 0;
        for (let i = 0; i < chr2nd.length; i++) {
          if (i < markNstart[2] || i > markNend[2]) {
            tempchr2nd.push(chr2nd[i]);
            tempchrorig.push(chrorig[i]);
            tempchrorient.push(chrorient[i]);
            tempdataphase.push(dataphase[i]);
            tempmarkerpos.push(markerpos[i]);
          } else if (flag1 === 0) {
            flag1 = 1;
            for (let j = chr2nd.length - 1; j >= 0; j--) {
              if (j >= markNstart[2] && j <= markNend[2]) {
                tempchr2nd.push(chr2nd[j]);
                tempchrorig.push(chrorig[j]);
                if (chrorient[j] === "+") {
                  tempchrorient.push("-");
                } else if (chrorient[j] === "-") {
                  tempchrorient.push("+");
                } else {
                  tempchrorient.push("na");
                }
                tempdataphase.push(dataphase[j]);
                tempmarkerpos.push(markerpos[j]);
              }
            }
          }
        }
        chr2nd = tempchr2nd;
        chrorig = tempchrorig;
        chrorient = tempchrorient;
        dataphase = tempdataphase;
        //console.log(dataphase);
        distphase = calcdist(dataphase);
        markerpos = tempmarkerpos;
        savedata();

        g.clear();
        for (let j = 0; j < distphase.length; j++) {
          for (let i = 0; i < distphase.length; i++) {
            g.beginFill(16777215 - (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)) - 256 * (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)));
            g.drawRect(j, i, 1, 1);
            g.endFill();
          }
        }

        spinner.classList.add('loaded');
      }, 1000)
    }
  }

  if (e.keyCode === 77) {
    console.log("m");
    if (markM === 0) {
      markM = 1;
      markMstart1[0] = -1;
      let tempReplace3 = chrorig[tempX1];
      for (let i = 0; i < chrorig.length; i++) {
        if (chrorig[i] === tempReplace3) {
          if (markMstart1[0] === -1) { markMstart1[0] = i };
          markMend1[0] = i;
        }
      }
    } else if (markM === 1) {
      markM = 2;
      markMstart1[1] = -1;
      let tempReplace3 = chrorig[tempX1];
      for (let i = 0; i < chrorig.length; i++) {
        if (chrorig[i] === tempReplace3) {
          if (markMstart1[1] === -1) { markMstart1[1] = i };
          markMend1[1] = i;
        }
      }
      if (markMstart1[0] < markMstart1[1]) {
        markMstart1[2] = markMstart1[0];
        markMend1[2] = markMend1[1];
      } else {
        markMstart1[2] = markMstart1[1];
        markMend1[2] = markMend1[0];
      }
      console.log("M1:" + markMstart1 + ", " + markMend1);
    } else if (markM === 2) {
      markM = 3;
      markMstart2[0] = -1;
      let tempReplace3 = chrorig[tempX1];
      for (let i = 0; i < chrorig.length; i++) {
        if (chrorig[i] === tempReplace3) {
          if (markMstart2[0] === -1) { markMstart2[0] = i };
          markMend2[0] = i;
        }
      }
    } else if (markM === 3) {
      const spinner = document.getElementById('loading');
      spinner.classList.remove('loaded');
      setTimeout(function () {
        markM = 0;
        markMstart2[1] = -1;
        let tempReplace3 = chrorig[tempX1];
        for (let i = 0; i < chrorig.length; i++) {
          if (chrorig[i] === tempReplace3) {
            if (markMstart2[1] === -1) { markMstart2[1] = i };
            markMend2[1] = i;
          }
        }
        if (markMstart2[0] < markMstart2[1]) {
          markMstart2[2] = markMstart2[0];
          markMend2[2] = markMend2[1];
        } else {
          markMstart2[2] = markMstart2[1];
          markMend2[2] = markMend2[0];
        }
        console.log("M2:" + markMstart2 + ", " + markMend2);

        g2.lineStyle(1, 0xFF0000);
        g2.drawRect(markMstart1[2], markMstart1[2], markMend1[2] - markMstart1[2] + 1, markMend1[2] - markMstart1[2] + 1);
        g2.drawRect(markMstart2[2], markMstart2[2], markMend2[2] - markMstart2[2] + 1, markMend2[2] - markMstart2[2] + 1);


        if (!(markMstart1[2] >= markMstart2[2] && markMstart1[2] <= markMend2[2])
          && !(markMend1[2] >= markMstart2[2] && markMend1[2] <= markMend2[2])
          && !(markMstart2[2] >= markMstart1[2] && markMstart2[2] <= markMend1[2])
          && !(markMend2[2] >= markMstart1[2] && markMend2[2] <= markMend1[2])) {
          let tempchr2nd = [];
          let tempchrorig = [];
          let tempchrorient = [];
          let tempdataphase = [];
          let tempmarkerpos = [];
          let flag1 = 0;
          let flag2 = 0;
          for (let i = 0; i < chr2nd.length; i++) {
            if (!(i >= markMstart1[2] && i <= markMend1[2]) && !(i >= markMstart2[2] && i <= markMend2[2])) {
              tempchr2nd.push(chr2nd[i]);
              tempchrorig.push(chrorig[i]);
              tempchrorient.push(chrorient[i]);
              tempdataphase.push(dataphase[i]);
              tempmarkerpos.push(markerpos[i]);
            } else if (i >= markMstart2[2] && i <= markMend2[2] && flag2 === 0) {
              flag2 = 1;
              for (let j = 0; j < chr2nd.length; j++) {
                if (j >= markMstart1[2] && j <= markMend1[2]) {
                  tempchr2nd.push(chr2nd[j]);
                  tempchrorig.push(chrorig[j]);
                  tempchrorient.push(chrorient[j]);
                  tempdataphase.push(dataphase[j]);
                  tempmarkerpos.push(markerpos[j]);
                }
              }
            } else if (i >= markMstart1[2] && i <= markMend1[2] && flag1 === 0) {
              flag1 = 1;
              for (let j = 0; j < chr2nd.length; j++) {
                if (j >= markMstart2[2] && j <= markMend2[2]) {
                  tempchr2nd.push(chr2nd[j]);
                  tempchrorig.push(chrorig[j]);
                  tempchrorient.push(chrorient[j]);
                  tempdataphase.push(dataphase[j]);
                  tempmarkerpos.push(markerpos[j]);
                }
              }
            }
          }
          chr2nd = tempchr2nd;
          chrorig = tempchrorig;
          chrorient = tempchrorient;
          dataphase = tempdataphase;
          //console.log(dataphase);
          distphase = calcdist(dataphase);
          markerpos = tempmarkerpos;
          savedata();
          g.clear();
          for (let j = 0; j < distphase.length; j++) {
            for (let i = 0; i < distphase.length; i++) {
              g.beginFill(16777215 - (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)) - 256 * (parseInt((parseFloat(distphase[j][i]) - 0.5) * 2 * 255)));
              g.drawRect(j, i, 1, 1);
              g.endFill();
            }
          }

        }

        spinner.classList.add('loaded');
      }, 1000);
    }
  }

  if (e.keyCode === 69) {
    console.log("e");
    if (markN === 0) {
      markN = 1;
      markNstart[0] = -1;
      let tempReplace3 = chrorig[tempX1];
      for (let i = 0; i < chrorig.length; i++) {
        if (chrorig[i] === tempReplace3) {
          if (markNstart[0] === -1) { markNstart[0] = i };
          markNend[0] = i;
        }
      }
    } else if (markN === 1) {
      markN = 0;
      markNstart[1] = -1;
      let tempReplace3 = chrorig[tempX1];
      for (let i = 0; i < chrorig.length; i++) {
        if (chrorig[i] === tempReplace3) {
          if (markNstart[1] === -1) { markNstart[1] = i };
          markNend[1] = i;
        }
      }
      if (markNstart[0] < markNstart[1]) {
        markNstart[2] = markNstart[0];
        markNend[2] = markNend[1];
      } else {
        markNstart[2] = markNstart[1];
        markNend[2] = markNend[0];
      }
      console.log(markNstart + ", " + markNend);
      g2.lineStyle(1, 0xFF0000);
      g2.drawRect(markNstart, markNstart, markNend - markNstart + 1, markNend - markNstart + 1);

      newchrname = window.prompt("新しい染色体名を入力してください", "");
      for (let i = 0; i < chr2nd.length; i++) {
        if (i >= markNstart[2] && i <= markNend[2]) {
          chr2nd[i] = newchrname;
        }
      }
      savedata();
    }
  }
});

$(pixijscanvas).on({
  'mouseenter': function (element) {
    $("body").append('<div class="tooltip" id="mytooltip1"></div>');
    $(".tooltip").css({ top: element.pageY + 20, left: element.pageX + 10 });
  },
  'mousemove': function (element) {
    $(".tooltip").css({ top: element.pageY + 20, left: element.pageX + 10 });
  },
  'mouseleave': function () {
    $(".tooltip").remove();
  }
});


function savedata() {
  let str = "";
  for (let i = 0; i < chr2nd.length; i++) {
    str += chr2nd[i] + "\t" + chr2nd[i];
    if (chrorient[i] === "+" || chrorient[i] === "-") {
      str += "\t+\t" + chrorient[i];
    } else {
      str += "\tna\t" + chrorient[i];
    }
    str += "\t" + chrorig[i] + "\t" + markerpos[i];
    for (let j = 0; j < dataphase[i].length; j++) {
      if (dataphase[i][j] === 1) {
        str += "\t1";
      } else if (dataphase[i][j] === -1) {
        str += "\t0";
      } else {
        str += "\t-1";
      }
    }
    str += "\n";
  }
  var blob = new Blob([str], { type: "text/csv" }); //配列に上記の文字列(str)を設定
  let link = document.querySelector('#savedata');
  link.href = URL.createObjectURL(blob);
  link.download = "savedate.txt";
  //link.click();
}
