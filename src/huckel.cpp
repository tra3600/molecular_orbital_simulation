#include "huckel.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace mo {

namespace {
HuckelSystem chain(const std::string& name, const std::string& desc, int n) {
    HuckelSystem s;
    s.name = name;
    s.description = desc;
    for (int i = 0; i < n; ++i) {
        HuckelAtom a;
        a.x = i * 0.87;
        a.y = (i % 2) * 0.5;
        s.atoms.push_back(a);
        if (i) s.bonds.push_back({i - 1, i, 1.0});
    }
    return s;
}
HuckelSystem ring(const std::string& name, const std::string& desc, int n) {
    HuckelSystem s;
    s.name = name;
    s.description = desc;
    double R = 0.5 / std::sin(M_PI / n);
    for (int i = 0; i < n; ++i) {
        HuckelAtom a;
        double t = M_PI / 2 + 2 * M_PI * i / n;
        a.x = R * std::cos(t);
        a.y = R * std::sin(t);
        s.atoms.push_back(a);
        s.bonds.push_back({i, (i + 1) % n, 1.0});
    }
    return s;
}
}  // namespace

std::vector<HuckelSystem> huckelLibrary() {
    std::vector<HuckelSystem> v;
    v.push_back(chain("ethylene", "Éthène : 2 centres, 2 électrons π", 2));
    v.push_back(chain("allyl", "Allyle (radical neutre ; --charge ±1 pour cation/anion)", 3));
    v.back().atoms[0].electrons = 1;
    v.push_back(chain("butadiene", "Buta-1,3-diène : conjugaison linéaire", 4));
    v.push_back(chain("hexatriene", "Hexa-1,3,5-triène", 6));
    v.push_back(ring("cyclobutadiene", "Cyclobutadiène : antiaromatique (4n électrons π)", 4));
    v.push_back(ring("cyclopentadienyl", "Radical cyclopentadiényle (--charge -1 : anion aromatique)", 5));
    v.push_back(ring("benzene", "Benzène : aromatique (4n+2 électrons π)", 6));
    v.push_back(ring("cyclooctatetraene", "Cyclooctatétraène (modèle plan)", 8));
    {
        HuckelSystem s = ring("pyridine", "Pyridine : N pyridinique (h=0,5 ; k=1)", 6);
        s.atoms[0].sym = "N";
        s.atoms[0].h = 0.5;
        v.push_back(s);
    }
    {
        HuckelSystem s = ring("pyrrole", "Pyrrole : N pyrrolique à 2 électrons π (h=1,5 ; k=0,8)", 5);
        s.atoms[0].sym = "N";
        s.atoms[0].h = 1.5;
        s.atoms[0].electrons = 2;
        s.bonds[0].k = 0.8;
        s.bonds[4].k = 0.8;
        v.push_back(s);
    }
    {
        HuckelSystem s = chain("acroleine", "Acroléine CH2=CH–CH=O (O : h=1 ; k=1)", 4);
        s.atoms[3].sym = "O";
        s.atoms[3].h = 1.0;
        v.push_back(s);
    }
    {
        HuckelSystem s;
        s.name = "naphthalene";
        s.description = "Naphtalène : 10 centres, deux cycles fusionnés";
        const double c = 0.866;
        double xy[10][2] = {{0, .5}, {c, 1}, {2 * c, .5}, {2 * c, -.5}, {c, -1}, {0, -.5}, {-c, -1}, {-2 * c, -.5}, {-2 * c, .5}, {-c, 1}};
        for (auto& p : xy) { HuckelAtom a; a.x = p[0]; a.y = p[1]; s.atoms.push_back(a); }
        for (int i = 0; i < 10; ++i) s.bonds.push_back({i, (i + 1) % 10, 1.0});
        s.bonds.push_back({0, 5, 1.0});
        v.push_back(s);
    }
    return v;
}

bool findHuckel(const std::string& name, HuckelSystem& out) {
    for (auto& s : huckelLibrary())
        if (s.name == name) { out = s; return true; }
    return false;
}

HuckelResult solveHuckel(const HuckelSystem& s) {
    int n = int(s.atoms.size());
    Mat H(n, n);
    for (int i = 0; i < n; ++i) H(i, i) = -s.atoms[i].h;  // on diagonalise −(h, k) : valeurs propres croissantes = x décroissants
    for (auto& b : s.bonds) H(b.a, b.b) = H(b.b, b.a) = -b.k;
    std::vector<double> w;
    Mat V;
    jacobiEigen(H, w, V);
    HuckelResult r;
    r.x.resize(n);
    for (int i = 0; i < n; ++i) r.x[i] = -w[i];
    r.C = V;
    int ne = -s.charge;
    for (auto& a : s.atoms) ne += a.electrons;
    r.occ.assign(n, 0.0);
    int left = ne;
    // remplissage aufbau ; les niveaux dégénérés partiels sont répartis équitablement
    for (int i = 0; i < n && left > 0;) {
        int j = i;
        while (j < n && std::fabs(r.x[j] - r.x[i]) < 1e-6) ++j;
        int g = j - i, cap = 2 * g, put = std::min(left, cap);
        for (int k = i; k < j; ++k) r.occ[k] = double(put) / g;
        left -= put;
        i = j;
    }
    r.charges.assign(n, 0.0);
    r.bondOrders = Mat(n, n);
    for (int i = 0; i < n; ++i) {
        double pop = 0;
        for (int k = 0; k < n; ++k) pop += r.occ[k] * r.C(i, k) * r.C(i, k);
        r.charges[i] = s.atoms[i].electrons - pop;
        for (int j = 0; j < n; ++j) {
            double p = 0;
            for (int k = 0; k < n; ++k) p += r.occ[k] * r.C(i, k) * r.C(j, k);
            r.bondOrders(i, j) = p;
        }
    }
    for (int k = 0; k < n; ++k) r.piEnergy += r.occ[k] * r.x[k];
    bool allC = true;
    for (auto& a : s.atoms) if (a.sym != "C") allC = false;
    if (allC) {
        r.hasDelocalization = true;
        r.delocalization = r.piEnergy - ne * 1.0;  // référence : ne électrons dans des liaisons éthylène (x = 1)
    }
    return r;
}

void printHuckel(std::ostream& os, const HuckelSystem& s, const HuckelResult& r) {
    int n = int(s.atoms.size());
    os << "\n=== Hückel : " << s.name << " ===\n" << s.description << "\n";
    os << "Énergies E = α + xβ (β < 0 : x grand = plus stable)\n";
    os << std::fixed << std::setprecision(4);
    for (int k = n - 1; k >= 0; --k) {
        os << "  ψ" << std::setw(2) << k + 1 << "  x = " << std::setw(8) << r.x[k] << "   occ " << std::setprecision(1) << r.occ[k] << std::setprecision(4) << "   [";
        for (int i = 0; i < n; ++i) os << std::setw(8) << r.C(i, k);
        os << " ]\n";
    }
    os << "Énergie π = " << n * 0 + r.piEnergy << " β";
    if (r.hasDelocalization) os << "   énergie de délocalisation = " << r.delocalization << " β";
    os << "\nCharges π : ";
    for (int i = 0; i < n; ++i) os << s.atoms[i].sym << i + 1 << "=" << std::showpos << r.charges[i] << std::noshowpos << " ";
    os << "\nOrdres de liaison π : ";
    for (auto& b : s.bonds) os << b.a + 1 << "–" << b.b + 1 << "=" << r.bondOrders(b.a, b.b) << " ";
    os << "\n" << std::defaultfloat;
}

static const char* HUCKEL_TEMPLATE = R"HTML(<!doctype html><html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Hückel</title>
<style>
:root{--bg:#0f1420;--card:#171e2e;--line:#2a3550;--fg:#e6ebf5;--mut:#8b97b3;--acc:#7cc4ff}
body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,sans-serif;padding:18px 24px}
h1{margin:0 0 4px}h2{font-size:15px;color:var(--acc);margin:0 0 8px}.sub{color:var(--mut);margin-bottom:12px}
.wrap{display:grid;grid-template-columns:minmax(260px,330px) 1fr;gap:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px}
.mos{display:grid;grid-template-columns:repeat(auto-fill,minmax(170px,1fr));gap:10px}
.mo{background:#0c111b;border:1px solid var(--line);border-radius:10px;padding:6px;text-align:center;cursor:pointer}
.mo.sel{border-color:#ffd34d}.mo small{color:var(--mut)}
.lvl{cursor:pointer}.lvl.sel line{stroke:#ffd34d!important}
.kpis{display:flex;gap:10px;flex-wrap:wrap;margin-bottom:12px}.kpi{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:6px 12px}.kpi b{display:block;font-size:16px}.kpi span{color:var(--mut);font-size:12px}
table{border-collapse:collapse}td,th{padding:2px 10px;border-bottom:1px solid var(--line);text-align:right}th{color:var(--mut);font-weight:500}
</style></head><body>
<h1 id="t"></h1><div class="sub" id="d"></div><div class="kpis" id="k"></div>
<div class="wrap"><div class="card"><h2>Niveaux π (E = α + xβ)</h2><div id="diag"></div></div>
<div class="card"><h2>Orbitales moléculaires π — cercle rouge : signe +, bleu : signe −, rayon ∝ |c|</h2><div class="mos" id="mos"></div>
<h2 style="margin-top:14px">Charges et ordres de liaison π</h2><div id="tabs"></div></div></div>
<script>
const D=__DATA__;const $=i=>document.getElementById(i);let sel=D.x.length-1;
$("t").textContent="Méthode de Hückel — "+D.name;$("d").textContent=D.desc;
$("k").innerHTML=[["Énergie π",D.pi.toFixed(4)+" β"],D.deloc===null?null:["Délocalisation",D.deloc.toFixed(4)+" β"],["Électrons π",D.ne]].filter(Boolean).map(k=>`<div class="kpi"><b>${k[1]}</b><span>${k[0]}</span></div>`).join("");
const n=D.x.length;
function mol(k,size){let minx=1e9,maxx=-1e9,miny=1e9,maxy=-1e9;D.atoms.forEach(a=>{minx=Math.min(minx,a.x);maxx=Math.max(maxx,a.x);miny=Math.min(miny,a.y);maxy=Math.max(maxy,a.y)});
 const sc=size/Math.max(maxx-minx+1.6,maxy-miny+1.6,2.6),cx=(minx+maxx)/2,cy=(miny+maxy)/2;
 const P=a=>[size/2+(a.x-cx)*sc,size/2-(a.y-cy)*sc];let s=`<svg viewBox="0 0 ${size} ${size}" width="100%">`;
 D.bonds.forEach(b=>{const p=P(D.atoms[b[0]]),q=P(D.atoms[b[1]]);s+=`<line x1="${p[0]}" y1="${p[1]}" x2="${q[0]}" y2="${q[1]}" stroke="#56627f" stroke-width="3"/>`});
 D.atoms.forEach((a,i)=>{const p=P(a),c=D.C[k][i],r=Math.abs(c)*sc*0.62*1.6;
  s+=`<circle cx="${p[0]}" cy="${p[1]}" r="${Math.max(0.5,r)}" fill="${c>=0?"#ff5a4d":"#4f8dff"}" fill-opacity=".75" stroke="#fff" stroke-opacity=".5"/>`;
  if(a.s!="C")s+=`<text x="${p[0]}" y="${p[1]+4}" text-anchor="middle" fill="#fff" font-size="11" font-weight="bold">${a.s}</text>`});
 return s+"</svg>";}
function draw(){
 let s=`<svg viewBox="0 0 300 ${60+n*34}" width="100%">`;const lo=Math.min(...D.x),hi=Math.max(...D.x),H=n*34+20;
 const y=x=>20+(hi-x)/Math.max(1e-6,hi-lo)*H;
 D.x.forEach((x,k)=>{let dup=0,tot=0;D.x.forEach((q,j)=>{if(Math.abs(q-x)<1e-6){if(j<k)dup++;tot++}});
  const x0=150+(dup-(tot-1)/2)*52-22,yy=y(x),occ=D.occ[k];
  const ar=occ>=2?"↑↓":occ>=1?"↑":occ>0?"½":"";
  s+=`<g class="lvl${k==sel?" sel":""}" data-k="${k}"><rect x="${x0-3}" y="${yy-10}" width="50" height="18" fill="transparent"/><line x1="${x0}" x2="${x0+44}" y1="${yy}" y2="${yy}" stroke="${occ>0?"#ff7a59":"#5aa0ff"}" stroke-width="3" stroke-linecap="round"/>
  <text x="${x0+22}" y="${yy-4}" text-anchor="middle" fill="#e6ebf5" font-size="12">${ar}</text>
  ${dup==0?`<text x="${x0-8}" y="${yy+4}" text-anchor="end" fill="#8b97b3" font-size="11">${(x>=0?"α+":"α−")+Math.abs(x).toFixed(3)}β</text>`:""}
  ${dup==tot-1?`<text x="${x0+50}" y="${yy+4}" fill="#8b97b3" font-size="11">${D.x.map((q,j)=>Math.abs(q-x)<1e-6?"ψ"+(j+1):null).filter(Boolean).join(" ")}</text>`:""}</g>`});
 $("diag").innerHTML=s+"</svg>";$("diag").querySelectorAll(".lvl").forEach(g=>g.onclick=()=>{sel=+g.dataset.k;draw();mos()});}
function mos(){$("mos").innerHTML=D.x.map((x,k)=>k).reverse().map(k=>`<div class="mo${k==sel?" sel":""}" data-k="${k}">${mol(k,170)}<small>ψ${k+1} · x=${D.x[k].toFixed(3)} · occ ${D.occ[k]}</small></div>`).join("");
 $("mos").querySelectorAll(".mo").forEach(g=>g.onclick=()=>{sel=+g.dataset.k;draw();mos()});}
$("tabs").innerHTML="<table><tr><th>atome</th><th>charge π</th></tr>"+D.atoms.map((a,i)=>`<tr><td>${a.s}${i+1}</td><td>${D.q[i]>=0?"+":"−"}${Math.abs(D.q[i]).toFixed(3)}</td></tr>`).join("")+"</table><br><table><tr><th>liaison</th><th>ordre π</th></tr>"+D.bonds.map(b=>`<tr><td>${b[0]+1}–${b[1]+1}</td><td>${b[2].toFixed(3)}</td></tr>`).join("")+"</table>";
draw();mos();
</script></body></html>
)HTML";

std::string huckelHTML(const HuckelSystem& s, const HuckelResult& r) {
    int n = int(s.atoms.size());
    auto num = [](double v) { std::ostringstream o; o << std::setprecision(8) << v; return o.str(); };
    std::ostringstream j;
    int ne = -s.charge;
    for (auto& a : s.atoms) ne += a.electrons;
    j << "{\"name\":\"" << s.name << "\",\"desc\":\"" << s.description << "\",\"ne\":" << ne << ",\"pi\":" << num(r.piEnergy)
      << ",\"deloc\":" << (r.hasDelocalization ? num(r.delocalization) : "null") << ",\"x\":[";
    for (int k = 0; k < n; ++k) j << (k ? "," : "") << num(r.x[k]);
    j << "],\"occ\":[";
    for (int k = 0; k < n; ++k) j << (k ? "," : "") << num(r.occ[k]);
    j << "],\"C\":[";
    for (int k = 0; k < n; ++k) {
        j << (k ? "," : "") << "[";
        // orientation de signe : plus grand coefficient positif
        double big = 0;
        for (int i = 0; i < n; ++i) if (std::fabs(r.C(i, k)) > std::fabs(big) + 1e-9) big = r.C(i, k);
        double sg = big < 0 ? -1 : 1;
        for (int i = 0; i < n; ++i) j << (i ? "," : "") << num(sg * r.C(i, k));
        j << "]";
    }
    j << "],\"atoms\":[";
    for (int i = 0; i < n; ++i) j << (i ? "," : "") << "{\"s\":\"" << s.atoms[i].sym << "\",\"x\":" << s.atoms[i].x << ",\"y\":" << s.atoms[i].y << "}";
    j << "],\"bonds\":[";
    for (size_t b = 0; b < s.bonds.size(); ++b) j << (b ? "," : "") << "[" << s.bonds[b].a << "," << s.bonds[b].b << "," << num(r.bondOrders(s.bonds[b].a, s.bonds[b].b)) << "]";
    j << "],\"q\":[";
    for (int i = 0; i < n; ++i) j << (i ? "," : "") << num(r.charges[i]);
    j << "]}";
    std::string html = HUCKEL_TEMPLATE;
    html.replace(html.find("__DATA__"), 8, j.str());
    return html;
}

}  // namespace mo
