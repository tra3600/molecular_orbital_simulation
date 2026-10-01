// Page HTML autonome (aucune dépendance externe) du rapport interactif.
#pragma once

namespace mo {
static const char* HTML_REPORT_TEMPLATE = R"HTML(<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Orbitales moléculaires</title>
<style>
:root{--bg:#0f1420;--card:#171e2e;--line:#2a3550;--fg:#e6ebf5;--mut:#8b97b3;--pos:#ff5a4d;--neg:#4f8dff;--acc:#7cc4ff}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,"Segoe UI",Roboto,sans-serif}
header{padding:18px 24px 6px}
h1{margin:0;font-size:24px;letter-spacing:.2px}
h2{margin:0 0 8px;font-size:15px;color:var(--acc);font-weight:600}
.sub{color:var(--mut);margin-top:4px}
.kpis{display:flex;flex-wrap:wrap;gap:10px;padding:10px 24px}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:8px 14px;min-width:130px}
.kpi b{display:block;font-size:17px}.kpi span{color:var(--mut);font-size:12px}
main{display:grid;grid-template-columns:minmax(270px,340px) minmax(380px,1fr) minmax(280px,380px);gap:14px;padding:6px 24px 14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px;min-width:0}
canvas{width:100%;height:auto;border-radius:8px;background:radial-gradient(circle at 50% 40%,#1c2640,#0c111b);touch-action:none;cursor:grab;display:block}
.controls{display:flex;flex-wrap:wrap;gap:10px 16px;margin-top:10px;align-items:center;color:var(--mut);font-size:13px}
.controls input[type=range]{width:110px;vertical-align:middle}
button{background:#22304f;color:var(--fg);border:1px solid var(--line);border-radius:7px;padding:4px 10px;cursor:pointer}
button:hover{background:#2b3d66}
.lvl{cursor:pointer}.lvl:hover line{stroke-width:5}.lvl.sel line{stroke:#ffd34d!important;stroke-width:5}
.lvl text{fill:var(--mut);font-size:11px}.lvl .arrows{fill:var(--fg);font-size:13px}
.legend span{display:inline-block;width:11px;height:11px;border-radius:50%;margin:0 4px -1px 8px}
.slices{display:grid;grid-template-columns:1fr 1fr;gap:8px}
.slices figure{margin:0}.slices img{width:100%;border-radius:6px;border:1px solid var(--line);image-rendering:auto}
figcaption{font-size:11px;color:var(--mut);text-align:center}
.bar{display:flex;align-items:center;gap:8px;margin:2px 0;font-size:12px}
.bar i{display:block;height:9px;border-radius:3px}.bar .l{width:76px;color:var(--mut)}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:14px;padding:0 24px 30px}
table{border-collapse:collapse;width:100%}td,th{padding:3px 8px;border-bottom:1px solid var(--line);text-align:right}
th:first-child,td:first-child{text-align:left}th{color:var(--mut);font-weight:500}
.tag{background:#ffd34d;color:#222;border-radius:5px;padding:0 6px;font-size:11px;margin-left:6px}
footer{padding:0 24px 24px;color:var(--mut);font-size:12px}
@media(max-width:1000px){main{grid-template-columns:1fr}}
</style></head><body>
<header><h1 id="title"></h1><div class="sub" id="sub"></div></header>
<div class="kpis" id="kpis"></div>
<main>
 <section class="card"><h2>Niveaux d'énergie <small style="color:var(--mut);font-weight:400">(cliquez sur un niveau)</small></h2><div id="diag"></div></section>
 <section class="card"><h2 id="orbTitle"></h2>
  <canvas id="cv" width="760" height="560"></canvas>
  <div class="controls">
   <span class="legend"><span style="background:var(--pos)"></span>lobe +<span style="background:var(--neg)"></span>lobe −</span>
   <label><input type="checkbox" id="auto" checked> rotation</label>
   <label>taille <input type="range" id="psz" min="1" max="8" step="0.5" value="3.5"></label>
   <label>opacité <input type="range" id="alp" min="0.1" max="1" step="0.05" value="0.55"></label>
   <label><input type="checkbox" id="showdens"> densité totale</label>
   <button id="reset">réinitialiser</button>
  </div>
  <div class="sub" style="margin-top:6px">Le nuage de points échantillonne |ψ|² : glisser = tourner, molette = zoom, ←/→ = orbitale précédente/suivante.</div>
 </section>
 <section class="card"><h2>Orbitale sélectionnée</h2><div id="orbInfo"></div><div class="slices" id="slices"></div></section>
</main>
<div class="grid2">
 <section class="card"><h2>Atomes</h2><table id="atomTab"></table></section>
 <section class="card"><h2>Liaisons (indice de Mayer)</h2><table id="bondTab"></table></section>
 <section class="card"><h2>Convergence SCF</h2><canvas id="conv" width="460" height="220" style="cursor:default"></canvas></section>
 <section class="card"><h2>Densité électronique totale</h2><div class="slices" id="dens"></div></section>
</div>
<footer>Rapport généré par molecular_orbital_simulation — Hartree–Fock / STO-3G, orbitales moléculaires par combinaison linéaire d'orbitales atomiques (LCAO). Les fichiers .cube du dossier s'ouvrent dans Avogadro, VMD, Jmol ou IboView.</footer>
<script>
const D = __DATA__;
const CPK={H:"#ffffff",He:"#d9ffff",Li:"#cc80ff",Be:"#c2ff00",B:"#ffb5b5",C:"#8a8a8a",N:"#3050f8",O:"#ff2a1a",F:"#90e050",Ne:"#b3e3f5"};
const $=id=>document.getElementById(id);
const EV=27.211386246;
let sel=0,yaw=0.6,pitch=-0.35,zoom=1,auto=true;

// ---------- en-tête
$("title").textContent=D.title+"  ·  "+D.formula;
$("sub").textContent=(D.desc?D.desc+" — ":"")+D.method+", charge "+D.charge+", multiplicité "+D.mult+", "+D.nbf+" fonctions de base";
const kp=[["Énergie totale",D.energy.toFixed(6)+" Eh"],["Énergie nucléaire",D.eNuc.toFixed(4)+" Eh"],
 ["Gap HOMO–LUMO",D.gap.toFixed(2)+" eV"],["Moment dipolaire",D.charge==0?D.dipole.toFixed(3)+" D":"—"],
 ["SCF",(D.converged?"convergé":"NON convergé")+" ("+D.iterations+" it.)"]];
if(D.uhf)kp.push(["⟨S²⟩",D.s2.toFixed(3)]);
$("kpis").innerHTML=kp.map(k=>`<div class="kpi"><b>${k[1]}</b><span>${k[0]}</span></div>`).join("");

// ---------- diagramme d'énergie
const f=e=>e>=-1.5?e:-1.5-0.17*Math.log(1+(-1.5-e));
function drawDiagram(){
  const top=34,bot=30,GAP=17;
  let lo=1e9,hi=-1e9;D.orbitals.forEach(o=>{lo=Math.min(lo,f(o.eps));hi=Math.max(hi,f(o.eps));});
  if(hi-lo<1e-6)hi=lo+1;
  const spins=D.uhf?[1,2]:[0];
  const W=D.uhf?420:320;
  // groupes de niveaux dégénérés par spin, espacés d'au moins GAP pixels
  const cols=spins.map(sp=>{
    const list=D.orbitals.map((o,i)=>({o,i})).filter(p=>p.o.spin==sp).sort((a,b)=>b.o.eps-a.o.eps);
    const groups=[];list.forEach(p=>{const g=groups[groups.length-1];if(g&&Math.abs(g[0].o.eps-p.o.eps)<2e-3)g.push(p);else groups.push([p]);});
    let prev=-1e9;groups.forEach(g=>{const y0=top+(hi-f(g[0].o.eps))/(hi-lo)*480;g.y=Math.max(y0,prev+GAP);prev=g.y;});
    return {sp,groups,end:prev};
  });
  const H=Math.max(520,Math.max(...cols.map(c=>c.end))+bot);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  cols.forEach(c=>{
    const cx=D.uhf?(c.sp==1?W*0.27:W*0.73):W*0.45;
    if(D.uhf)s+=`<text x="${cx}" y="18" text-anchor="middle" fill="#8b97b3" font-size="12">${c.sp==1?"α ↑":"β ↓"}</text>`;
    c.groups.forEach(g=>{
      const tot=g.length,sw=Math.min(46,(W-120)/tot),xs=cx-(tot*sw)/2+3;
      g.forEach((p,k)=>{const o=p.o,x0=xs+k*sw,yy=g.y,col=o.occ>0?"#ff7a59":"#5aa0ff";
        const ar=o.occ>=2?"↑↓":(o.occ==1?(c.sp==2?"↓":"↑"):"");
        s+=`<g class="lvl${p.i==sel?" sel":""}" data-i="${p.i}"><rect x="${x0-3}" y="${yy-10}" width="${sw}" height="18" fill="transparent"/>
        <line x1="${x0}" x2="${x0+sw-8}" y1="${yy}" y2="${yy}" stroke="${col}" stroke-width="3" stroke-linecap="round"/>
        <text class="arrows" x="${x0+(sw-8)/2}" y="${yy-3}" text-anchor="middle">${ar}</text></g>`;});
      s+=`<text x="${xs-6}" y="${g.y+4}" text-anchor="end" fill="#8b97b3" font-size="11">${g[0].o.eps.toFixed(2)}</text>`;
      const tg=g.find(p=>p.o.tag);
      s+=`<text x="${xs+tot*sw}" y="${g.y+4}" fill="#8b97b3" font-size="11">${tot>3?g[0].o.name+"–"+g[tot-1].o.name:g.map(p=>p.o.name).join(" ")}${tg?" · "+tg.o.tag:""}</text>`;
    });
  });
  s+=`<text x="6" y="${H-8}" fill="#8b97b3" font-size="10">énergie (Eh) — axe comprimé sous −1,5</text></svg>`;
  $("diag").innerHTML=s;
  $("diag").querySelectorAll(".lvl").forEach(g=>g.onclick=()=>select(+g.dataset.i));
}

// ---------- sélection
function select(i){
  sel=(i+D.orbitals.length)%D.orbitals.length;
  const o=D.orbitals[sel];
  $("orbTitle").innerHTML=`Orbitale ${o.name} <small style="color:var(--mut);font-weight:400">E = ${o.eps.toFixed(4)} Eh = ${(o.eps*EV).toFixed(2)} eV</small>`+(o.tag?`<span class="tag">${o.tag}</span>`:"");
  const mx=Math.max(...o.comp.map(c=>Math.abs(c[1])),1e-9);
  $("orbInfo").innerHTML="<div style='color:var(--mut);margin-bottom:4px'>Combinaison linéaire d'orbitales atomiques :</div>"+
   o.comp.map(c=>`<div class="bar"><span class="l">${c[0]}</span><i style="width:${Math.abs(c[1])/mx*120}px;background:${c[1]>=0?"var(--pos)":"var(--neg)"}"></i><span>${c[1]>=0?"+":"−"}${Math.abs(c[1]).toFixed(3)}</span></div>`).join("");
  $("slices").innerHTML=`<figure><img src="data:image/png;base64,${o.sa}"><figcaption>coupe 1 (axe principal × axe 2)</figcaption></figure><figure><img src="data:image/png;base64,${o.sb}"><figcaption>coupe 2 (axe principal × axe 3)</figcaption></figure>`;
  drawDiagram();
}

// ---------- vue 3D
const cv=$("cv"),ctx=cv.getContext("2d");
let cen=[0,0,0],rmax=3;
(function(){const n=D.atoms.length;D.atoms.forEach(a=>{cen[0]+=a.x/n;cen[1]+=a.y/n;cen[2]+=a.z/n;});
 rmax=1;D.atoms.forEach(a=>{rmax=Math.max(rmax,Math.hypot(a.x-cen[0],a.y-cen[1],a.z-cen[2]));});})();
function rot(x,y,z){const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);
  let x1=cy*x+sy*z,z1=-sy*x+cy*z;let y1=cp*y-sp*z1,z2=sp*y+cp*z1;return [x1,y1,z2];}
function render(){
  const W=cv.width,H=cv.height;ctx.clearRect(0,0,W,H);
  const scale=zoom*0.40*Math.min(W,H)/(rmax+2.6);
  const o=D.orbitals[sel],items=[];
  const alpha=+$("alp").value,ps=+$("psz").value;
  if($("showdens").checked){/* densité : toutes les orbitales occupées */
    D.orbitals.forEach(q=>{if(q.occ>0&&(q.spin!=2||D.uhf))add(q,Math.min(1,q.occ)*1,true);});
  } else add(o,1,false);
  function add(q,w,dens){const p=q.pts,step=dens?Math.max(1,Math.round(p.length/4/700))*4:4;
    for(let i=0;i<p.length;i+=step){const r=rot(p[i]/100-cen[0],p[i+1]/100-cen[1],p[i+2]/100-cen[2]);
      items.push({z:r[2],x:W/2+r[0]*scale,y:H/2-r[1]*scale,k:0,s:dens?0:p[i+3]});}}
  D.atoms.forEach((a,i)=>{const r=rot(a.x-cen[0],a.y-cen[1],a.z-cen[2]);items.push({z:r[2],x:W/2+r[0]*scale,y:H/2-r[1]*scale,k:1,a:a,i:i,rr:r});});
  D.bonds.forEach(b=>{if(b[2]<0.35)return;const A=D.atoms[b[0]],B=D.atoms[b[1]];
    const r1=rot(A.x-cen[0],A.y-cen[1],A.z-cen[2]),r2=rot(B.x-cen[0],B.y-cen[1],B.z-cen[2]);
    items.push({z:(r1[2]+r2[2])/2-0.01,k:2,x1:W/2+r1[0]*scale,y1:H/2-r1[1]*scale,x2:W/2+r2[0]*scale,y2:H/2-r2[1]*scale,o:b[2]});});
  items.sort((a,b)=>a.z-b.z);
  const rad=ps*scale/40;
  items.forEach(it=>{
    if(it.k==0){ctx.globalAlpha=alpha*(0.55+0.45*(it.z/(rmax+4)+0.5>1?1:Math.max(0,it.z/(rmax+4)+0.5)));
      ctx.fillStyle=it.s==0?"#ffd34d":(it.s>0?"#ff5a4d":"#4f8dff");ctx.beginPath();ctx.arc(it.x,it.y,rad,0,6.2832);ctx.fill();}
    else if(it.k==2){ctx.globalAlpha=0.9;ctx.strokeStyle="#9aa6c4";ctx.lineWidth=Math.min(9,2+2.2*it.o)*zoom;ctx.lineCap="round";
      ctx.beginPath();ctx.moveTo(it.x1,it.y1);ctx.lineTo(it.x2,it.y2);ctx.stroke();}
    else{ctx.globalAlpha=0.95;const R=Math.max(4,it.a.r*0.42*scale);
      const g=ctx.createRadialGradient(it.x-R*.3,it.y-R*.3,R*.1,it.x,it.y,R);const c=CPK[it.a.s]||"#bbb";
      g.addColorStop(0,"#ffffff");g.addColorStop(.35,c);g.addColorStop(1,"#10141f");
      ctx.fillStyle=g;ctx.beginPath();ctx.arc(it.x,it.y,R,0,6.2832);ctx.fill();
      ctx.globalAlpha=1;ctx.fillStyle="#e6ebf5";ctx.font="12px system-ui";ctx.fillText(it.a.s+(it.i+1),it.x+R+3,it.y-R*0.6);}
  });
  ctx.globalAlpha=1;
}
let drag=null;
cv.addEventListener("pointerdown",e=>{drag=[e.clientX,e.clientY];cv.setPointerCapture(e.pointerId);cv.style.cursor="grabbing";});
cv.addEventListener("pointermove",e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*0.008;pitch=Math.max(-1.55,Math.min(1.55,pitch+(e.clientY-drag[1])*0.008));drag=[e.clientX,e.clientY];});
cv.addEventListener("pointerup",()=>{drag=null;cv.style.cursor="grab";});
cv.addEventListener("wheel",e=>{e.preventDefault();zoom=Math.max(0.3,Math.min(5,zoom*(e.deltaY<0?1.1:0.91)));},{passive:false});
$("reset").onclick=()=>{yaw=0.6;pitch=-0.35;zoom=1;};
$("auto").onchange=e=>auto=e.target.checked;
document.addEventListener("keydown",e=>{if(e.key=="ArrowRight")select(sel+1);if(e.key=="ArrowLeft")select(sel-1);});
function loop(){if(auto&&!drag)yaw+=0.006;render();requestAnimationFrame(loop);}

// ---------- tableaux & convergence
$("atomTab").innerHTML="<tr><th>atome</th><th>charge de Mulliken</th>"+(D.uhf?"<th>spin</th>":"")+"</tr>"+
 D.atoms.map((a,i)=>`<tr><td>${a.s}${i+1}</td><td>${a.q>=0?"+":"−"}${Math.abs(a.q).toFixed(3)}</td>${D.uhf?`<td>${a.spin.toFixed(3)}</td>`:""}</tr>`).join("");
$("bondTab").innerHTML="<tr><th>liaison</th><th>longueur (Å)</th><th>indice</th></tr>"+
 (D.bonds.map(b=>`<tr><td>${D.atoms[b[0]].s}${b[0]+1}–${D.atoms[b[1]].s}${b[1]+1}</td><td>${b[3].toFixed(3)}</td><td>${b[2].toFixed(3)}</td></tr>`).join("")||"<tr><td colspan=3>aucune</td></tr>");
$("dens").innerHTML=`<figure><img src="data:image/png;base64,${D.density[0]}"><figcaption>coupe 1</figcaption></figure><figure><img src="data:image/png;base64,${D.density[1]}"><figcaption>coupe 2</figcaption></figure>`;
(function(){const c=$("conv"),x=c.getContext("2d"),h=D.history,n=h.length,Ef=h[n-1];
 x.fillStyle="#0c111b";x.fillRect(0,0,c.width,c.height);
 const v=h.map(e=>Math.log10(Math.max(1e-12,Math.abs(e-Ef))));const lo=-12,hi=Math.max(1,...v);
 x.strokeStyle="#2a3550";x.fillStyle="#8b97b3";x.font="11px system-ui";
 for(let k=-12;k<=hi;k+=3){const yy=10+(hi-k)/(hi-lo)*(c.height-40);x.beginPath();x.moveTo(40,yy);x.lineTo(c.width-10,yy);x.stroke();x.fillText("1e"+k,2,yy+4);}
 x.strokeStyle="#7cc4ff";x.lineWidth=2;x.beginPath();
 v.forEach((q,i)=>{const px=40+i/Math.max(1,n-1)*(c.width-55),py=10+(hi-q)/(hi-lo)*(c.height-40);i?x.lineTo(px,py):x.moveTo(px,py);});x.stroke();
 x.fillText("itération → (|E − E_final| en hartree)",50,c.height-8);})();
select(0);loop();
</script></body></html>
)HTML";
}  // namespace mo
